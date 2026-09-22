#include <md_gisaxs.h>

#include <core/md_fft.h>
#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define GISAXS_MAX_CLASSES 4
#define GISAXS_DEFAULT_MAX_SLICES 1024
#define GISAXS_DEFAULT_OVERSAMPLING 2.0

// ------------------------------------------------------------------------------------------------
// Small helpers
// ------------------------------------------------------------------------------------------------

typedef struct cplx_t {
    double re, im;
} cplx_t;

static inline cplx_t c_make(double re, double im) { cplx_t c = {re, im}; return c; }
static inline cplx_t c_add(cplx_t a, cplx_t b) { return c_make(a.re + b.re, a.im + b.im); }
static inline cplx_t c_sub(cplx_t a, cplx_t b) { return c_make(a.re - b.re, a.im - b.im); }
static inline cplx_t c_mul(cplx_t a, cplx_t b) { return c_make(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re); }
static inline cplx_t c_scale(cplx_t a, double s) { return c_make(a.re * s, a.im * s); }
static inline cplx_t c_div(cplx_t a, cplx_t b) {
    const double d = b.re * b.re + b.im * b.im;
    return c_make((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
}
// Principal square root (Re >= 0, Im has the sign of the input imaginary part, Im >= 0 for negative reals)
static inline cplx_t c_sqrt(cplx_t a) {
    const double r = hypot(a.re, a.im);
    const double re = sqrt(MAX(0.0, (r + a.re) * 0.5));
    double im = sqrt(MAX(0.0, (r - a.re) * 0.5));
    if (a.im < 0.0) im = -im;
    return c_make(re, im);
}
static inline cplx_t c_exp(cplx_t a) {
    const double e = exp(a.re);
    return c_make(e * cos(a.im), e * sin(a.im));
}
static inline cplx_t c_sin(cplx_t a) {
    // sin(x + iy) = sin x cosh y + i cos x sinh y
    return c_make(sin(a.re) * cosh(a.im), cos(a.re) * sinh(a.im));
}

static inline double sinc(double x) {
    return fabs(x) < 1.0e-6 ? 1.0 - x * x / 6.0 : sin(x) / x;
}

static inline cplx_t c_sinc(cplx_t x) {
    if (hypot(x.re, x.im) < 1.0e-6) {
        cplx_t x2 = c_mul(x, x);
        return c_make(1.0 - x2.re / 6.0, -x2.im / 6.0);
    }
    return c_div(c_sin(x), x);
}

// Cubic B-spline weights for the 4 nodes i-1, i, i+1, i+2 where i = floor(u), t = u - i
static inline void bspline4(float t, float w[4]) {
    const float t2 = t * t;
    const float t3 = t2 * t;
    const float s = 1.0f - t;
    w[0] = s * s * s * (1.0f / 6.0f);
    w[1] = (3.0f * t3 - 6.0f * t2 + 4.0f) * (1.0f / 6.0f);
    w[2] = (-3.0f * t3 + 3.0f * t2 + 3.0f * t + 1.0f) * (1.0f / 6.0f);
    w[3] = t3 * (1.0f / 6.0f);
}

// Fourier transform of the cubic B-spline assignment window for a phase x = q * h / 2
static inline double bspline4_window(double x) {
    const double s = sinc(x);
    return (s * s) * (s * s);
}

static inline cplx_t c_bspline4_window(cplx_t x) {
    const cplx_t s = c_sinc(x);
    const cplx_t s2 = c_mul(s, s);
    return c_mul(s2, s2);
}

static inline size_t packed_row_offset(size_t row, size_t n) {
    // Offset of row 'row' in a row major packed upper triangular matrix of dimension n
    return row * n - (row * (row - 1)) / 2;
}

static void* aligned_alloc_zero(size_t bytes) {
    void* ptr = md_fft_alloc(bytes);
    if (ptr) MEMSET(ptr, 0, bytes);
    return ptr;
}

// ------------------------------------------------------------------------------------------------
// Context
// ------------------------------------------------------------------------------------------------

typedef struct particle_t {
    float u;        // z in slice units: (z - z0) / dz
    float x, y;     // In-plane position in grid units [0, nx), [0, ny)
    float w;        // Weight
    uint32_t cls;   // Class index
} particle_t;

struct md_gisaxs_t {
    struct md_allocator_i* alloc;

    // Particles sorted by u
    size_t      num_particles;
    particle_t* particles;
    float*      particle_u;     // Sorted copy of u for binary search
    double      particle_z_min;
    double      particle_z_max;

    size_t num_classes;
    double class_sigma[GISAXS_MAX_CLASSES];

    double box_x, box_y, area;
    int    nx, ny, nkx;
    double dx, dy;
    int    kx_count;
    md_fft_2d_t* fft;

    size_t  num_slices;
    double  dz, z0;
    double* slice_z;
    double* profile;

    // Extended dimension: num_classes * num_slices
    size_t dim;
    size_t packed_size;

    size_t    num_rings;
    double    dq;
    double*   ring_q;
    unsigned* ring_count;
    uint32_t* ring_offset;      // num_rings + 1, offsets into points

    size_t    num_points;
    uint32_t* point_index;      // index into r2c output (ky * nkx + kx)
    int16_t*  point_k;          // (kx, ky_signed) pairs, used by the reference implementation
    float*    point_scale;      // num_points * num_classes: in-plane Gaussian / B-spline window

    float* spectra;             // num_points * 2 * dim   [p][re(dim) im(dim)]
    float* matrix;              // num_rings * packed_size

    size_t spectra_bytes;
    size_t matrix_bytes;
};

static int cmp_particle(const void* a, const void* b) {
    const float ua = ((const particle_t*)a)->u;
    const float ub = ((const particle_t*)b)->u;
    return (ua > ub) - (ua < ub);
}

static int cmp_float(const void* a, const void* b) {
    const float fa = *(const float*)a;
    const float fb = *(const float*)b;
    return (fa > fb) - (fa < fb);
}

// Assign particles to at most GISAXS_MAX_CLASSES classes of equal sigma
static size_t build_classes(double out_sigma[GISAXS_MAX_CLASSES], uint32_t* out_cls, const md_gisaxs_input_t* in, const uint8_t* include) {
    const size_t N = in->count;
    if (!in->sigma) {
        out_sigma[0] = MAX(0.0, (double)in->sigma_uniform);
        for (size_t i = 0; i < N; ++i) out_cls[i] = 0;
        return 1;
    }

    // Collect unique values (quantized to 1e-4 Å)
    float uniq[GISAXS_MAX_CLASSES + 1];
    size_t num_uniq = 0;
    bool overflow = false;
    for (size_t i = 0; i < N && !overflow; ++i) {
        if (!include[i]) continue;
        const float s = roundf(MAX(0.0f, in->sigma[i]) * 1.0e4f) * 1.0e-4f;
        size_t j = 0;
        for (; j < num_uniq; ++j) {
            if (uniq[j] == s) break;
        }
        if (j == num_uniq) {
            if (num_uniq == GISAXS_MAX_CLASSES) {
                overflow = true;
            } else {
                uniq[num_uniq++] = s;
            }
        }
    }

    if (!overflow) {
        if (num_uniq == 0) {
            uniq[0] = 0.0f;
            num_uniq = 1;
        }
        for (size_t j = 0; j < num_uniq; ++j) out_sigma[j] = uniq[j];
        for (size_t i = 0; i < N; ++i) {
            out_cls[i] = 0;
            if (!include[i]) continue;
            const float s = roundf(MAX(0.0f, in->sigma[i]) * 1.0e4f) * 1.0e-4f;
            for (uint32_t j = 0; j < num_uniq; ++j) {
                if (uniq[j] == s) { out_cls[i] = j; break; }
            }
        }
        return num_uniq;
    }

    // Too many distinct widths: quantize into classes of equal population (by sorted sigma),
    // class sigma is the RMS sigma of its members.
    MD_LOG_INFO("GISAXS: more than %i distinct Gaussian widths, grouping into %i classes", GISAXS_MAX_CLASSES, GISAXS_MAX_CLASSES);
    size_t num_incl = 0;
    for (size_t i = 0; i < N; ++i) num_incl += include[i] ? 1 : 0;
    float* sorted = (float*)malloc(sizeof(float) * MAX(num_incl, 1));
    size_t n = 0;
    for (size_t i = 0; i < N; ++i) if (include[i]) sorted[n++] = MAX(0.0f, in->sigma[i]);
    qsort(sorted, n, sizeof(float), cmp_float);
    float bounds[GISAXS_MAX_CLASSES];
    for (int c = 0; c < GISAXS_MAX_CLASSES; ++c) {
        size_t idx = MIN(n - 1, ((size_t)(c + 1) * n) / GISAXS_MAX_CLASSES);
        bounds[c] = sorted[idx];
    }
    bounds[GISAXS_MAX_CLASSES - 1] = FLT_MAX;
    free(sorted);

    double sum2[GISAXS_MAX_CLASSES] = {0};
    size_t cnt[GISAXS_MAX_CLASSES] = {0};
    for (size_t i = 0; i < N; ++i) {
        out_cls[i] = 0;
        if (!include[i]) continue;
        const float s = MAX(0.0f, in->sigma[i]);
        uint32_t c = 0;
        while (c < GISAXS_MAX_CLASSES - 1 && s > bounds[c]) ++c;
        out_cls[i] = c;
        sum2[c] += (double)s * s;
        cnt[c] += 1;
    }
    // Compact empty classes
    uint32_t remap[GISAXS_MAX_CLASSES];
    size_t num = 0;
    for (int c = 0; c < GISAXS_MAX_CLASSES; ++c) {
        if (cnt[c]) {
            remap[c] = (uint32_t)num;
            out_sigma[num] = sqrt(sum2[c] / (double)cnt[c]);
            ++num;
        } else {
            remap[c] = 0;
        }
    }
    for (size_t i = 0; i < N; ++i) out_cls[i] = remap[out_cls[i]];
    return MAX(num, 1);
}

md_gisaxs_t* md_gisaxs_create(const md_gisaxs_input_t* in, const md_gisaxs_params_t* params, struct md_allocator_i* alloc) {
    if (!in || !params || !alloc) {
        MD_LOG_ERROR("GISAXS: invalid arguments");
        return NULL;
    }
    if (!in->x || !in->y || !in->z) {
        MD_LOG_ERROR("GISAXS: missing coordinates");
        return NULL;
    }
    if (!(in->box_x > 0.0) || !(in->box_y > 0.0)) {
        MD_LOG_ERROR("GISAXS: a periodic XY box is required");
        return NULL;
    }
    if (!(params->q_par_max > 0.0) || !(params->q_z_max > 0.0)) {
        MD_LOG_ERROR("GISAXS: q_par_max and q_z_max must be positive");
        return NULL;
    }

    const size_t N = in->count;
    const double os = params->oversampling > 0.0 ? MAX(params->oversampling, 1.0) : GISAXS_DEFAULT_OVERSAMPLING;
    const size_t max_slices = params->max_slices ? params->max_slices : GISAXS_DEFAULT_MAX_SLICES;
    const bool   use_range = params->z_max > params->z_min;

    md_gisaxs_t* ctx = (md_gisaxs_t*)md_alloc(alloc, sizeof(md_gisaxs_t));
    MEMSET(ctx, 0, sizeof(md_gisaxs_t));
    ctx->alloc = alloc;
    ctx->box_x = in->box_x;
    ctx->box_y = in->box_y;
    ctx->area  = in->box_x * in->box_y;

    // --- Particle selection and z extent ---
    uint8_t* include = (uint8_t*)malloc(MAX(N, 1));
    double zmin = DBL_MAX, zmax = -DBL_MAX;
    size_t num_incl = 0;
    size_t num_clipped = 0;
    for (size_t i = 0; i < N; ++i) {
        const float w = in->weight ? in->weight[i] : 1.0f;
        const double z = in->z[i];
        bool inc = (w != 0.0f) && isfinite(z) && isfinite(in->x[i]) && isfinite(in->y[i]);
        if (inc && use_range && (z < params->z_min || z > params->z_max)) {
            inc = false;
            num_clipped += 1;
        }
        include[i] = inc;
        if (inc) {
            zmin = MIN(zmin, z);
            zmax = MAX(zmax, z);
            num_incl += 1;
        }
    }
    if (num_clipped) {
        MD_LOG_INFO("GISAXS: %zu particles outside the z range were excluded", num_clipped);
    }
    if (num_incl == 0) {
        MD_LOG_ERROR("GISAXS: no particles to process");
        free(include);
        md_gisaxs_destroy(ctx);
        return NULL;
    }
    ctx->particle_z_min = zmin;
    ctx->particle_z_max = zmax;

    // --- Classes ---
    uint32_t* cls = (uint32_t*)malloc(sizeof(uint32_t) * MAX(N, 1));
    ctx->num_classes = build_classes(ctx->class_sigma, cls, in, include);

    // --- Slices ---
    double dz = params->dz > 0.0 ? params->dz : PI / (os * params->q_z_max);
    // B-spline support: a particle at u touches floor(u)-1 .. floor(u)+2
    double z0 = zmin - 2.0 * dz;
    size_t num_slices = (size_t)ceil((zmax - z0) / dz) + 3;
    if (num_slices > max_slices) {
        const double new_dz = (zmax - zmin) / (double)(max_slices - 5);
        MD_LOG_INFO("GISAXS: slice spacing increased from %.3f to %.3f Å to respect the maximum number of slices (%zu). Accuracy at high q_z is reduced.", dz, new_dz, max_slices);
        dz = new_dz;
        z0 = zmin - 2.0 * dz;
        num_slices = (size_t)ceil((zmax - z0) / dz) + 3;
    }
    ctx->dz = dz;
    ctx->z0 = z0;
    ctx->num_slices = num_slices;
    ctx->dim = ctx->num_classes * num_slices;
    ctx->packed_size = ctx->dim * (ctx->dim + 1) / 2;

    ctx->slice_z = (double*)md_alloc(alloc, sizeof(double) * num_slices);
    ctx->profile = (double*)md_alloc(alloc, sizeof(double) * num_slices);
    for (size_t k = 0; k < num_slices; ++k) {
        ctx->slice_z[k] = z0 + (double)k * dz;
        ctx->profile[k] = 0.0;
    }

    // --- In-plane grid ---
    const double q_max = params->q_par_max;
    ctx->nx = md_fft_valid_size((int)ceil(os * q_max * in->box_x / PI), true);
    ctx->ny = md_fft_valid_size((int)ceil(os * q_max * in->box_y / PI), true);
    ctx->nkx = ctx->nx / 2 + 1;
    ctx->dx = in->box_x / ctx->nx;
    ctx->dy = in->box_y / ctx->ny;
    ctx->fft = md_fft_2d_create(ctx->nx, ctx->ny);
    if (!ctx->fft) {
        MD_LOG_ERROR("GISAXS: failed to create FFT plan (%i x %i)", ctx->nx, ctx->ny);
        free(include);
        free(cls);
        md_gisaxs_destroy(ctx);
        return NULL;
    }

    const double dqx = 2.0 * PI / in->box_x;
    const double dqy = 2.0 * PI / in->box_y;
    const int kx_max = (int)floor(q_max / dqx);
    const int ky_max = (int)floor(q_max / dqy);
    ctx->kx_count = MIN(kx_max + 1, ctx->nkx);

    // --- Rings ---
    ctx->dq = MAX(dqx, dqy);
    ctx->num_rings = (size_t)floor(q_max / ctx->dq + 0.5);
    if (ctx->num_rings == 0) {
        MD_LOG_ERROR("GISAXS: q_par_max is smaller than the reciprocal grid spacing of the box");
        free(include);
        free(cls);
        md_gisaxs_destroy(ctx);
        return NULL;
    }

    // Half plane points: kx > 0 (all ky), kx == 0 (ky > 0). Each represents itself and its conjugate mirror.
    // Count per ring
    const size_t R = ctx->num_rings;
    ctx->ring_offset = (uint32_t*)md_alloc(alloc, sizeof(uint32_t) * (R + 1));
    ctx->ring_count  = (unsigned*)md_alloc(alloc, sizeof(unsigned) * R);
    ctx->ring_q      = (double*)md_alloc(alloc, sizeof(double) * R);
    MEMSET(ctx->ring_offset, 0, sizeof(uint32_t) * (R + 1));
    MEMSET(ctx->ring_count, 0, sizeof(unsigned) * R);
    MEMSET(ctx->ring_q, 0, sizeof(double) * R);

#define RING_OF(kx, ky, out_q, out_r) do { \
        const double _qx = (kx) * dqx, _qy = (ky) * dqy; \
        out_q = sqrt(_qx * _qx + _qy * _qy); \
        long _r = (long)floor(out_q / ctx->dq + 0.5); \
        if (_r < 1) _r = 1; \
        out_r = _r - 1; \
    } while(0)

    size_t num_points = 0;
    for (int kx = 0; kx <= kx_max && kx < ctx->nkx; ++kx) {
        for (int ky = -ky_max; ky <= ky_max; ++ky) {
            if (kx == 0 && ky <= 0) continue;
            if (ky <= -(ctx->ny / 2) || ky >= ctx->ny / 2) continue;
            double q; long r;
            RING_OF(kx, ky, q, r);
            if (q > q_max || r >= (long)R) continue;
            ctx->ring_offset[r + 1] += 1;
            ctx->ring_q[r] += q;
            num_points += 1;
        }
    }
    for (size_t r = 0; r < R; ++r) {
        ctx->ring_count[r] = 2 * ctx->ring_offset[r + 1];
        ctx->ring_q[r] = ctx->ring_offset[r + 1] ? ctx->ring_q[r] / ctx->ring_offset[r + 1] : (r + 1) * ctx->dq;
        ctx->ring_offset[r + 1] += ctx->ring_offset[r];
    }
    ctx->num_points = num_points;

    const size_t C = ctx->num_classes;
    ctx->point_index = (uint32_t*)md_alloc(alloc, sizeof(uint32_t) * MAX(num_points, 1));
    ctx->point_k     = (int16_t*) md_alloc(alloc, sizeof(int16_t) * 2 * MAX(num_points, 1));
    ctx->point_scale = (float*)   md_alloc(alloc, sizeof(float) * C * MAX(num_points, 1));
    {
        uint32_t* fill = (uint32_t*)malloc(sizeof(uint32_t) * R);
        MEMCPY(fill, ctx->ring_offset, sizeof(uint32_t) * R);
        for (int kx = 0; kx <= kx_max && kx < ctx->nkx; ++kx) {
            for (int ky = -ky_max; ky <= ky_max; ++ky) {
                if (kx == 0 && ky <= 0) continue;
                if (ky <= -(ctx->ny / 2) || ky >= ctx->ny / 2) continue;
                double q; long r;
                RING_OF(kx, ky, q, r);
                if (q > q_max || r >= (long)R) continue;
                const uint32_t p = fill[r]++;
                const int ky_idx = ky < 0 ? ky + ctx->ny : ky;
                ctx->point_index[p] = (uint32_t)(ky_idx * ctx->nkx + kx);
                ctx->point_k[2 * p + 0] = (int16_t)kx;
                ctx->point_k[2 * p + 1] = (int16_t)ky;
                const double wx = bspline4_window(PI * kx / ctx->nx);
                const double wy = bspline4_window(PI * ky / ctx->ny);
                for (size_t c = 0; c < C; ++c) {
                    const double s = ctx->class_sigma[c];
                    ctx->point_scale[p * C + c] = (float)(exp(-0.5 * q * q * s * s) / (wx * wy));
                }
            }
        }
        free(fill);
    }
#undef RING_OF

    // --- Particles (sorted by u) ---
    ctx->particles  = (particle_t*)md_alloc(alloc, sizeof(particle_t) * num_incl);
    ctx->particle_u = (float*)md_alloc(alloc, sizeof(float) * num_incl);
    {
        size_t n = 0;
        const double inv_dx = 1.0 / ctx->dx;
        const double inv_dy = 1.0 / ctx->dy;
        const double inv_dz = 1.0 / dz;
        for (size_t i = 0; i < N; ++i) {
            if (!include[i]) continue;
            double x = fmod(in->x[i], in->box_x); if (x < 0) x += in->box_x;
            double y = fmod(in->y[i], in->box_y); if (y < 0) y += in->box_y;
            particle_t p;
            p.u = (float)((in->z[i] - z0) * inv_dz);
            p.x = (float)(x * inv_dx);
            p.y = (float)(y * inv_dy);
            if (p.x >= (float)ctx->nx) p.x -= (float)ctx->nx;
            if (p.y >= (float)ctx->ny) p.y -= (float)ctx->ny;
            p.w = in->weight ? in->weight[i] : 1.0f;
            p.cls = cls[i];
            ctx->particles[n++] = p;
        }
        ctx->num_particles = n;
        qsort(ctx->particles, n, sizeof(particle_t), cmp_particle);
        for (size_t i = 0; i < n; ++i) ctx->particle_u[i] = ctx->particles[i].u;
    }
    free(include);
    free(cls);

    // --- Storage ---
    ctx->spectra_bytes = sizeof(float) * 2 * ctx->dim * MAX(num_points, 1);
    ctx->matrix_bytes  = sizeof(float) * ctx->packed_size * R;
    ctx->spectra = (float*)aligned_alloc_zero(ctx->spectra_bytes);
    ctx->matrix  = (float*)aligned_alloc_zero(ctx->matrix_bytes);
    if (!ctx->spectra || !ctx->matrix) {
        MD_LOG_ERROR("GISAXS: failed to allocate %.1f MB for spectra and %.1f MB for ring matrices",
            ctx->spectra_bytes / (1024.0 * 1024.0), ctx->matrix_bytes / (1024.0 * 1024.0));
        md_gisaxs_destroy(ctx);
        return NULL;
    }

    return ctx;
}

void md_gisaxs_destroy(md_gisaxs_t* ctx) {
    if (!ctx) return;
    struct md_allocator_i* alloc = ctx->alloc;
    const size_t R = ctx->num_rings;
    const size_t P = MAX(ctx->num_points, 1);
    if (ctx->particles)   md_free(alloc, ctx->particles,   sizeof(particle_t) * ctx->num_particles);
    if (ctx->particle_u)  md_free(alloc, ctx->particle_u,  sizeof(float) * ctx->num_particles);
    if (ctx->slice_z)     md_free(alloc, ctx->slice_z,     sizeof(double) * ctx->num_slices);
    if (ctx->profile)     md_free(alloc, ctx->profile,     sizeof(double) * ctx->num_slices);
    if (ctx->ring_offset) md_free(alloc, ctx->ring_offset, sizeof(uint32_t) * (R + 1));
    if (ctx->ring_count)  md_free(alloc, ctx->ring_count,  sizeof(unsigned) * R);
    if (ctx->ring_q)      md_free(alloc, ctx->ring_q,      sizeof(double) * R);
    if (ctx->point_index) md_free(alloc, ctx->point_index, sizeof(uint32_t) * P);
    if (ctx->point_k)     md_free(alloc, ctx->point_k,     sizeof(int16_t) * 2 * P);
    if (ctx->point_scale) md_free(alloc, ctx->point_scale, sizeof(float) * ctx->num_classes * P);
    if (ctx->spectra)     md_fft_free(ctx->spectra);
    if (ctx->matrix)      md_fft_free(ctx->matrix);
    if (ctx->fft)         md_fft_2d_destroy(ctx->fft);
    md_free(alloc, ctx, sizeof(md_gisaxs_t));
}

void md_gisaxs_get_info(const md_gisaxs_t* ctx, md_gisaxs_info_t* info) {
    if (!info) return;
    MEMSET(info, 0, sizeof(md_gisaxs_info_t));
    if (!ctx) return;
    info->nx = ctx->nx;
    info->ny = ctx->ny;
    info->dx = ctx->dx;
    info->dy = ctx->dy;
    info->num_slices = ctx->num_slices;
    info->dz = ctx->dz;
    info->z0 = ctx->z0;
    info->num_rings = ctx->num_rings;
    info->dq_ring = ctx->dq;
    info->num_points = ctx->num_points;
    info->num_classes = ctx->num_classes;
    for (size_t c = 0; c < ctx->num_classes && c < 4; ++c) info->class_sigma[c] = ctx->class_sigma[c];
    info->num_particles = ctx->num_particles;
    info->area = ctx->area;
    info->spectra_bytes = ctx->spectra_bytes;
    info->matrix_bytes = ctx->matrix_bytes;
}

// ------------------------------------------------------------------------------------------------
// Slices
// ------------------------------------------------------------------------------------------------

static inline size_t grid_floats(const md_gisaxs_t* ctx) {
    return ALIGN_TO((size_t)ctx->nx * (size_t)ctx->ny, 16);
}

static inline size_t spec_floats(const md_gisaxs_t* ctx) {
    return ALIGN_TO((size_t)ctx->nkx * (size_t)ctx->ny * 2, 16);
}

size_t md_gisaxs_slice_scratch_bytes(const md_gisaxs_t* ctx) {
    if (!ctx) return 0;
    const size_t floats = ctx->num_classes * grid_floats(ctx) + spec_floats(ctx) + ALIGN_TO(md_fft_2d_scratch_size(ctx->fft), 16);
    return floats * sizeof(float);
}

static size_t lower_bound_f(const float* arr, size_t n, float v) {
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        if (arr[mid] < v) lo = mid + 1; else hi = mid;
    }
    return lo;
}

void md_gisaxs_compute_slices(md_gisaxs_t* ctx, size_t beg, size_t end, void* scratch) {
    ASSERT(ctx);
    ASSERT(scratch);
    end = MIN(end, ctx->num_slices);
    if (beg >= end) return;

    const int nx = ctx->nx;
    const int ny = ctx->ny;
    const size_t C = ctx->num_classes;
    const size_t S = ctx->num_slices;
    const size_t M = ctx->dim;
    const size_t P = ctx->num_points;

    float* grids = (float*)scratch;
    float* spec  = grids + C * grid_floats(ctx);
    float* fft_scratch = spec + spec_floats(ctx);

    for (size_t k = beg; k < end; ++k) {
        MEMSET(grids, 0, sizeof(float) * C * grid_floats(ctx));

        // Particles with u in (k - 2, k + 2] contribute to slice k (floor(u) in [k-2, k+1])
        const float u_lo = (float)k - 2.0f;
        const float u_hi = (float)k + 2.0f;
        size_t i = lower_bound_f(ctx->particle_u, ctx->num_particles, u_lo);
        double total = 0.0;
        bool   any = false;
        for (; i < ctx->num_particles; ++i) {
            const particle_t p = ctx->particles[i];
            if (p.u >= u_hi) break;
            const float fu = floorf(p.u);
            const int   iu = (int)fu;
            const int   node = (int)k - iu + 1;   // Which of the 4 spline nodes (iu-1 .. iu+2) slice k is
            if (node < 0 || node > 3) continue;
            float wz[4];
            bspline4(p.u - fu, wz);
            const float wk = p.w * wz[node];
            if (wk == 0.0f) continue;
            total += wk;
            any = true;

            float* grid = grids + p.cls * grid_floats(ctx);
            const float fx = floorf(p.x);
            const float fy = floorf(p.y);
            float wx[4], wy[4];
            bspline4(p.x - fx, wx);
            bspline4(p.y - fy, wy);
            int ix[4], iy[4];
            for (int j = 0; j < 4; ++j) {
                int a = (int)fx - 1 + j;
                int b = (int)fy - 1 + j;
                a = a < 0 ? a + nx : (a >= nx ? a - nx : a);
                b = b < 0 ? b + ny : (b >= ny ? b - ny : b);
                ix[j] = a;
                iy[j] = b;
            }
            for (int jy = 0; jy < 4; ++jy) {
                float* row = grid + (size_t)iy[jy] * nx;
                const float wyk = wy[jy] * wk;
                row[ix[0]] += wx[0] * wyk;
                row[ix[1]] += wx[1] * wyk;
                row[ix[2]] += wx[2] * wyk;
                row[ix[3]] += wx[3] * wyk;
            }
        }
        ctx->profile[k] = total / (ctx->area * ctx->dz);

        for (size_t c = 0; c < C; ++c) {
            const float* grid = grids + c * grid_floats(ctx);
            float* dst_base = ctx->spectra + c * S + k;
            if (!any) {
                for (size_t p = 0; p < P; ++p) {
                    dst_base[p * 2 * M + 0] = 0.0f;
                    dst_base[p * 2 * M + M] = 0.0f;
                }
                continue;
            }
            md_fft_2d_r2c(ctx->fft, grid, spec, ctx->kx_count, fft_scratch);
            for (size_t p = 0; p < P; ++p) {
                const uint32_t idx = ctx->point_index[p];
                const float s = ctx->point_scale[p * C + c];
                dst_base[p * 2 * M + 0] = spec[2 * idx + 0] * s;
                dst_base[p * 2 * M + M] = spec[2 * idx + 1] * s;
            }
        }
    }
}

// ------------------------------------------------------------------------------------------------
// Rings
// ------------------------------------------------------------------------------------------------

void md_gisaxs_compute_rings(md_gisaxs_t* ctx, size_t beg, size_t end) {
    ASSERT(ctx);
    end = MIN(end, ctx->num_rings);
    if (beg >= end) return;
    if (!ctx->spectra) {
        MD_LOG_ERROR("GISAXS: spectra have been released");
        return;
    }

    const size_t M = ctx->dim;
    double* acc = (double*)md_fft_alloc(sizeof(double) * ctx->packed_size);

    for (size_t r = beg; r < end; ++r) {
        MEMSET(acc, 0, sizeof(double) * ctx->packed_size);
        const uint32_t p_beg = ctx->ring_offset[r];
        const uint32_t p_end = ctx->ring_offset[r + 1];
        for (uint32_t p = p_beg; p < p_end; ++p) {
            const float* re = ctx->spectra + (size_t)p * 2 * M;
            const float* im = re + M;
            for (size_t m = 0; m < M; ++m) {
                const double rm = re[m];
                const double ii = im[m];
                if (rm == 0.0 && ii == 0.0) continue;
                double* row = acc + packed_row_offset(m, M) - m;
                for (size_t n = m; n < M; ++n) {
                    row[n] += rm * re[n] + ii * im[n];
                }
            }
        }
        const size_t count = p_end - p_beg;
        const double scl = count ? 1.0 / (double)count : 0.0;
        float* dst = ctx->matrix + r * ctx->packed_size;
        for (size_t i = 0; i < ctx->packed_size; ++i) {
            dst[i] = (float)(acc[i] * scl);
        }
    }

    md_fft_free(acc);
}

void md_gisaxs_release_spectra(md_gisaxs_t* ctx) {
    if (!ctx || !ctx->spectra) return;
    md_fft_free(ctx->spectra);
    ctx->spectra = NULL;
}

bool md_gisaxs_compute(md_gisaxs_t* ctx) {
    if (!ctx) return false;
    void* scratch = md_fft_alloc(md_gisaxs_slice_scratch_bytes(ctx));
    if (!scratch) return false;
    md_gisaxs_compute_slices(ctx, 0, ctx->num_slices, scratch);
    md_fft_free(scratch);
    md_gisaxs_compute_rings(ctx, 0, ctx->num_rings);
    md_gisaxs_release_spectra(ctx);
    return true;
}

// ------------------------------------------------------------------------------------------------
// Accessors
// ------------------------------------------------------------------------------------------------

size_t md_gisaxs_num_rings(const md_gisaxs_t* ctx)            { return ctx ? ctx->num_rings : 0; }
const double* md_gisaxs_ring_q(const md_gisaxs_t* ctx)         { return ctx ? ctx->ring_q : NULL; }
const unsigned* md_gisaxs_ring_count(const md_gisaxs_t* ctx)   { return ctx ? ctx->ring_count : NULL; }
size_t md_gisaxs_num_slices(const md_gisaxs_t* ctx)           { return ctx ? ctx->num_slices : 0; }
const double* md_gisaxs_slice_z(const md_gisaxs_t* ctx)        { return ctx ? ctx->slice_z : NULL; }
const double* md_gisaxs_slice_profile(const md_gisaxs_t* ctx)  { return ctx ? ctx->profile : NULL; }
double md_gisaxs_particle_z_min(const md_gisaxs_t* ctx)        { return ctx ? ctx->particle_z_min : 0.0; }
double md_gisaxs_particle_z_max(const md_gisaxs_t* ctx)        { return ctx ? ctx->particle_z_max : 0.0; }

// ------------------------------------------------------------------------------------------------
// Evaluation
// ------------------------------------------------------------------------------------------------

// Layered reference medium. Layers are indexed from the top: 0 = ambient (half space above h[0]),
// 1..N = film layers (graded), N+1 = substrate (half space below h[N]).
// h[j] is the bottom boundary of layer j (j = 0..N), d[j] the thickness of layer j (j = 1..N).
// The field in layer j is T_j exp(-i k_j (z - h_j)) + R_j exp(+i k_j (z - h_j)).
typedef struct stack_t {
    size_t  num_film;       // N
    cplx_t* sld;            // N+2, relative to the ambient (layer 0 is 0)
    double* h;              // N+1
    double* d;              // N+2 (d[0] and d[N+1] unused)
    double  roughness2;     // substrate roughness^2
    // Slice -> layer mapping
    uint32_t* slice_layer;  // num_slices
} stack_t;

typedef struct field_t {
    cplx_t* k;              // N+2
    cplx_t* T;              // N+1 (T in substrate not needed)
    cplx_t* R;              // N+1
} field_t;

static inline cplx_t c_rdiv_safe(cplx_t num, cplx_t den) {
    if (hypot(den.re, den.im) < 1.0e-300) return c_make(0.0, 0.0);
    return c_div(num, den);
}

static void stack_init(stack_t* st, const md_gisaxs_t* ctx, const md_gisaxs_model_t* m) {
    const size_t S = ctx->num_slices;
    const double dz = ctx->dz;
    size_t k0 = S;  // First slice (from the bottom) that belongs to the film
    if (m->graded) {
        for (size_t k = 0; k < S; ++k) {
            if (ctx->slice_z[k] + 0.5 * dz > m->z_substrate) { k0 = k; break; }
        }
    }
    const size_t N = S - k0;
    st->num_film = N;
    st->sld = (cplx_t*)md_fft_alloc(sizeof(cplx_t) * (N + 2));
    st->h   = (double*)md_fft_alloc(sizeof(double) * (N + 1));
    st->d   = (double*)md_fft_alloc(sizeof(double) * (N + 2));
    st->slice_layer = (uint32_t*)md_fft_alloc(sizeof(uint32_t) * MAX(S, 1));
    st->roughness2 = m->substrate_roughness * m->substrate_roughness;

    st->sld[0] = c_make(0.0, 0.0);
    st->d[0] = 0.0;
    if (N == 0) {
        st->h[0] = m->z_substrate;
    } else {
        st->h[0] = ctx->slice_z[S - 1] + 0.5 * dz;
        for (size_t j = 1; j <= N; ++j) {
            const size_t k = S - j;     // layer 1 = top slice
            const double top = st->h[j - 1];
            const double bot = (j == N) ? m->z_substrate : ctx->slice_z[k] - 0.5 * dz;
            st->h[j] = bot;
            st->d[j] = MAX(top - bot, 0.0);
            // SLD = re - i*abs, relative to the ambient
            st->sld[j] = c_make(m->profile_sld_scale * ctx->profile[k], -m->profile_abs_scale * ctx->profile[k]);
        }
    }
    st->sld[N + 1] = c_make(m->sld_substrate - m->sld_ambient, -m->sld_substrate_abs);
    st->d[N + 1] = 0.0;

    for (size_t k = 0; k < S; ++k) {
        if (N == 0) {
            st->slice_layer[k] = 0;
        } else if (k >= k0) {
            st->slice_layer[k] = (uint32_t)(S - k);
        } else {
            // Below the substrate interface: continue the field of the lowest film layer (only B-spline tails end up here)
            st->slice_layer[k] = (uint32_t)N;
        }
    }
}

static void stack_free(stack_t* st) {
    md_fft_free(st->sld);
    md_fft_free(st->h);
    md_fft_free(st->d);
    md_fft_free(st->slice_layer);
}

static void field_alloc(field_t* f, size_t N) {
    f->k = (cplx_t*)md_fft_alloc(sizeof(cplx_t) * (N + 2));
    f->T = (cplx_t*)md_fft_alloc(sizeof(cplx_t) * (N + 1));
    f->R = (cplx_t*)md_fft_alloc(sizeof(cplx_t) * (N + 1));
}

static void field_free(field_t* f) {
    md_fft_free(f->k);
    md_fft_free(f->T);
    md_fft_free(f->R);
}

// Parratt recursion for a wave with vertical wave number kz0 (> 0, ambient) incident from the top with unit amplitude.
// X (scratch, N+2) holds R/T at the bottom of each layer.
static void field_compute(field_t* f, cplx_t* X, const stack_t* st, double kz0) {
    const size_t N = st->num_film;
    // Vertical wave numbers: k_j^2 = kz0^2 - 4 pi dSLD_j
    for (size_t j = 0; j < N + 2; ++j) {
        const cplx_t s = st->sld[j];
        f->k[j] = c_sqrt(c_make(kz0 * kz0 - 4.0 * PI * s.re, -4.0 * PI * s.im));
    }
    f->k[0] = c_make(kz0, 0.0);

    // Upward recursion of the ratio X_j = R_j / T_j at the bottom of layer j
    X[N + 1] = c_make(0.0, 0.0);
    for (size_t jj = N + 1; jj-- > 0;) {
        const size_t j = jj;
        cplx_t r = c_rdiv_safe(c_sub(f->k[j], f->k[j + 1]), c_add(f->k[j], f->k[j + 1]));
        if (j == N && st->roughness2 > 0.0) {
            // Nevot-Croce
            r = c_mul(r, c_exp(c_scale(c_mul(f->k[j], f->k[j + 1]), -2.0 * st->roughness2)));
        }
        // Ratio of layer j+1 at its top boundary (z = h_j)
        cplx_t Xp = c_make(0.0, 0.0);
        if (j + 1 <= N) {
            const cplx_t ph = c_exp(c_mul(c_make(0.0, 2.0 * st->d[j + 1]), f->k[j + 1]));
            Xp = c_mul(X[j + 1], ph);
        }
        X[j] = c_rdiv_safe(c_add(r, Xp), c_add(c_make(1.0, 0.0), c_mul(r, Xp)));
    }

    // Downward amplitudes
    f->T[0] = c_make(1.0, 0.0);
    f->R[0] = X[0];
    for (size_t j = 0; j < N; ++j) {
        const cplx_t ph2 = c_exp(c_mul(c_make(0.0, 2.0 * st->d[j + 1]), f->k[j + 1]));
        const cplx_t Xp = c_mul(X[j + 1], ph2);
        // Continuity at h_j: T_j (1 + X_j) = T'_{j+1} (1 + X'_{j+1}), T_{j+1} = T'_{j+1} exp(i k d)
        const cplx_t Tp = c_rdiv_safe(c_mul(f->T[j], c_add(c_make(1.0, 0.0), X[j])), c_add(c_make(1.0, 0.0), Xp));
        const cplx_t ph = c_exp(c_mul(c_make(0.0, st->d[j + 1]), f->k[j + 1]));
        f->T[j + 1] = c_mul(Tp, ph);
        f->R[j + 1] = c_mul(X[j + 1], f->T[j + 1]);
    }
}

// exp(-Q^2 sigma^2 / 2) / W(Q dz / 2) for complex Q
static inline cplx_t z_factor(cplx_t Q, double sigma, double dz) {
    const cplx_t Q2 = c_mul(Q, Q);
    const cplx_t g = c_exp(c_scale(Q2, -0.5 * sigma * sigma));
    const cplx_t W = c_bspline4_window(c_scale(Q, 0.5 * dz));
    return c_rdiv_safe(g, W);
}

static void accumulate_intensity(const md_gisaxs_t* ctx, const double* cr, const double* ci, double scale, float* dst) {
    const size_t R = ctx->num_rings;
    const size_t M = ctx->dim;
    // I = sum_m S_mm |c_m|^2 + 2 sum_{m<n} S_mn (cr_m cr_n + ci_m ci_n)
    for (size_t r = 0; r < R; ++r) {
        const float* mat = ctx->matrix + r * ctx->packed_size;
        double I = 0.0;
        for (size_t m = 0; m < M; ++m) {
            const double a = cr[m];
            const double b = ci[m];
            if (a == 0.0 && b == 0.0) continue;
            const float* row = mat + packed_row_offset(m, M) - m;
            double sa = 0.0, sb = 0.0;
            for (size_t n = m + 1; n < M; ++n) {
                sa += row[n] * cr[n];
                sb += row[n] * ci[n];
            }
            I += row[m] * (a * a + b * b) + 2.0 * (a * sa + b * sb);
        }
        dst[r] = (float)(MAX(I, 0.0) * scale);
    }
}

void md_gisaxs_evaluate_range(const md_gisaxs_t* ctx, const md_gisaxs_model_t* model, const double* qz, size_t qz_beg, size_t qz_end, float* out) {
    ASSERT(ctx);
    ASSERT(model);
    ASSERT(qz);
    ASSERT(out);
    if (model->dwba && !(model->wavelength > 0.0)) {
        MD_LOG_ERROR("GISAXS: invalid wavelength");
        return;
    }

    const size_t R = ctx->num_rings;
    const size_t S = ctx->num_slices;
    const size_t C = ctx->num_classes;
    const size_t M = ctx->dim;
    const double dz = ctx->dz;
    const double scale = (model->intensity_scale != 0.0 ? model->intensity_scale : 1.0) / ctx->area;

    double* cr = (double*)md_fft_alloc(sizeof(double) * 2 * M);
    double* ci = cr + M;

    if (!model->dwba) {
        // Born approximation: a single term Q = q_z
        for (size_t iq = qz_beg; iq < qz_end; ++iq) {
            const cplx_t Q = c_make(qz[iq], 0.0);
            for (size_t c = 0; c < C; ++c) {
                const cplx_t f = z_factor(Q, ctx->class_sigma[c], dz);
                for (size_t k = 0; k < S; ++k) {
                    const double z = ctx->slice_z[k];
                    const cplx_t v = c_mul(f, c_make(cos(qz[iq] * z), -sin(qz[iq] * z)));
                    cr[c * S + k] = v.re;
                    ci[c * S + k] = v.im;
                }
            }
            accumulate_intensity(ctx, cr, ci, scale, out + iq * R);
        }
        md_fft_free(cr);
        return;
    }

    stack_t st;
    stack_init(&st, ctx, model);
    const size_t N = st.num_film;
    field_t fi, ff;
    field_alloc(&fi, N);
    field_alloc(&ff, N);
    cplx_t* X = (cplx_t*)md_fft_alloc(sizeof(cplx_t) * (N + 2));

    const double k0 = 2.0 * PI / model->wavelength;
    const double p = k0 * sin(model->alpha_i);   // incident, downward (ambient)
    field_compute(&fi, X, &st, p);

    for (size_t iq = qz_beg; iq < qz_end; ++iq) {
        float* dst = out + iq * R;
        const double q = qz[iq] - p;             // exit, upward (ambient)
        if (q < 0.0) {                           // Below the horizon
            MEMSET(dst, 0, sizeof(float) * R);
            continue;
        }
        field_compute(&ff, X, &st, q);

        for (size_t k = 0; k < S; ++k) {
            const uint32_t j = st.slice_layer[k];
            const cplx_t pj = fi.k[j];
            const cplx_t qj = ff.k[j];
            const double zl = ctx->slice_z[k] - st.h[j];
            cplx_t Q[4];
            cplx_t coef[4];
            Q[0] = c_add(pj, qj);                        coef[0] = c_mul(fi.T[j], ff.T[j]);
            Q[1] = c_sub(qj, pj);                        coef[1] = c_mul(fi.R[j], ff.T[j]);
            Q[2] = c_sub(pj, qj);                        coef[2] = c_mul(fi.T[j], ff.R[j]);
            Q[3] = c_scale(c_add(pj, qj), -1.0);         coef[3] = c_mul(fi.R[j], ff.R[j]);
            cplx_t tc[4];
            for (int t = 0; t < 4; ++t) {
                // coef * exp(-i Q zl), with complex Q = a + ib: exp(-i a zl + b zl)
                tc[t] = c_mul(coef[t], c_exp(c_make(Q[t].im * zl, -Q[t].re * zl)));
            }
            for (size_t c = 0; c < C; ++c) {
                cplx_t sum = c_make(0.0, 0.0);
                for (int t = 0; t < 4; ++t) {
                    sum = c_add(sum, c_mul(tc[t], z_factor(Q[t], ctx->class_sigma[c], dz)));
                }
                cr[c * S + k] = sum.re;
                ci[c * S + k] = sum.im;
            }
        }
        accumulate_intensity(ctx, cr, ci, scale, dst);
    }

    md_fft_free(X);
    field_free(&fi);
    field_free(&ff);
    stack_free(&st);
    md_fft_free(cr);
}

void md_gisaxs_evaluate(const md_gisaxs_t* ctx, const md_gisaxs_model_t* model, const double* qz, size_t num_qz, float* out) {
    md_gisaxs_evaluate_range(ctx, model, qz, 0, num_qz, out);
}

void md_gisaxs_reflectivity(const md_gisaxs_t* ctx, const md_gisaxs_model_t* model, const double* qz, size_t num_qz, double* out) {
    ASSERT(ctx);
    ASSERT(model);
    stack_t st;
    stack_init(&st, ctx, model);
    field_t f;
    field_alloc(&f, st.num_film);
    cplx_t* X = (cplx_t*)md_fft_alloc(sizeof(cplx_t) * (st.num_film + 2));
    for (size_t i = 0; i < num_qz; ++i) {
        const double kz = 0.5 * qz[i];
        if (kz <= 0.0) { out[i] = 1.0; continue; }
        field_compute(&f, X, &st, kz);
        out[i] = X[0].re * X[0].re + X[0].im * X[0].im;
    }
    md_fft_free(X);
    field_free(&f);
    stack_free(&st);
}

// ------------------------------------------------------------------------------------------------
// Reference
// ------------------------------------------------------------------------------------------------

void md_gisaxs_reference_born(const md_gisaxs_t* ctx, const md_gisaxs_input_t* in, const double* qz, size_t num_qz, double* out) {
    ASSERT(ctx);
    ASSERT(in);
    const size_t R = ctx->num_rings;
    const double dqx = 2.0 * PI / ctx->box_x;
    const double dqy = 2.0 * PI / ctx->box_y;

    for (size_t iq = 0; iq < num_qz; ++iq) {
        for (size_t r = 0; r < R; ++r) {
            double acc = 0.0;
            const uint32_t p_beg = ctx->ring_offset[r];
            const uint32_t p_end = ctx->ring_offset[r + 1];
            for (uint32_t p = p_beg; p < p_end; ++p) {
                const double qx = ctx->point_k[2 * p + 0] * dqx;
                const double qy = ctx->point_k[2 * p + 1] * dqy;
                const double q2 = qx * qx + qy * qy + qz[iq] * qz[iq];
                // Average over q_par and -q_par (full plane)
                for (int sgn = -1; sgn <= 1; sgn += 2) {
                    double re = 0.0, im = 0.0;
                    for (size_t j = 0; j < in->count; ++j) {
                        const double w = in->weight ? in->weight[j] : 1.0;
                        const double s = in->sigma ? in->sigma[j] : in->sigma_uniform;
                        const double amp = w * exp(-0.5 * q2 * s * s);
                        const double phase = -(sgn * (qx * in->x[j] + qy * in->y[j]) + qz[iq] * in->z[j]);
                        re += amp * cos(phase);
                        im += amp * sin(phase);
                    }
                    acc += 0.5 * (re * re + im * im);
                }
            }
            const size_t count = p_end - p_beg;
            out[iq * R + r] = count ? acc / (double)count / ctx->area : 0.0;
        }
    }
}

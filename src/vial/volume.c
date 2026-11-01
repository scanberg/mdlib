// Volumes: sdf() and density() (volume.h)

#include "volume.h"
#include "builtins.h"

#include <md_system.h>
#include <md_util.h>
#include <md_unitcell.h>
#include <core/md_spatial_acc.h>
#include <core/md_coord_stream.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_intrinsics.h>

#include <float.h>
#include <math.h>
#include <string.h>

#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#endif

bool vial_k_positions(vial_kctx_t* k, const vial_val_t* s, vial_val_t* out);

bool vial_volume_builtin(const char* name) {
    return name && (!strcmp(name, "sdf") || !strcmp(name, "density"));
}

double vial_volume_grid(const double len[3], double samples, uint32_t dims[3], bool* coarser) {
    bool c = false;
    double s = samples > 0.0 ? samples : VIAL_VOLUME_SAMPLES;
    for (int it = 0; it < 64; ++it) {
        uint64_t total = 1;
        for (int a = 0; a < 3; ++a) {
            // A length a whole number of steps long is that many voxels, not one more for the rounding of the product:
            // samples= written as 128 / 12 is a float
            double v = ceil(len[a] * s * (1.0 - 1.0e-6));
            if (!(v >= 1.0)) v = 1.0;
            if (v > 1.0e6) v = 1.0e6;
            const uint64_t d = ((uint64_t)v + 7) & ~(uint64_t)7;
            dims[a] = (uint32_t)d;
            total *= d;
        }
        if (total <= VIAL_VOLUME_MAX_VOXELS) break;
        s *= cbrt((double)VIAL_VOLUME_MAX_VOXELS / (double)total) * 0.999;
        c = true;
    }
    if (coarser) *coarser = c;
    return s;
}

// ### ARGUMENTS ###

typedef struct vargs_t {
    const vial_val_t* const* v;
    const str_t*      names;
    uint32_t          n;
} vargs_t;

// By parameter name when the check bound them, else the pos-th positional one; NULL when not given
static const vial_val_t* varg(const vargs_t* a, const char* name, uint32_t pos) {
    bool bound = false;
    for (uint32_t i = 0; i < a->n; ++i) {
        if (!a->names || !a->names[i].len) continue;
        if (str_eq(a->names[i], str_from_cstr(name))) return a->v[i];
        bound = true;
    }
    if (bound) return NULL;
    uint32_t p = 0;
    for (uint32_t i = 0; i < a->n; ++i) {
        if (a->names && a->names[i].len) continue;
        if (p++ == pos) return a->v[i];
    }
    return NULL;
}

static bool is_none(const vial_val_t* v) { return !v || v->dtype == VIAL_DT_NONE; }

static double num(const vial_val_t* v, size_t i) {
    switch (v->dtype) {
    case VIAL_DT_BOOL:    return ((const uint8_t*)v->data)[i];
    case VIAL_DT_INT:     return ((const int32_t*)v->data)[i];
    case VIAL_DT_FLOAT:   return ((const float*)v->data)[i];
    case VIAL_DT_FLOAT64: return ((const double*)v->data)[i];
    default:              return NAN;
    }
}

static bool is_num(const vial_val_t* v) {
    return v && v->data && (v->dtype == VIAL_DT_INT || v->dtype == VIAL_DT_FLOAT || v->dtype == VIAL_DT_FLOAT64 || v->dtype == VIAL_DT_BOOL);
}

// Three numbers: one for all three, a tuple of three, or float[3]
static bool three(const vial_val_t* v, double out[3]) {
    if (!v) return false;
    if (v->dtype == VIAL_DT_TUPLE) {
        if (v->num_elems != 3) return false;
        for (int a = 0; a < 3; ++a) {
            if (!is_num(&v->elems[a]) || vial_val_count(&v->elems[a]) != 1) return false;
            out[a] = num(&v->elems[a], 0);
        }
        return true;
    }
    if (!is_num(v)) return false;
    const size_t c = vial_val_count(v);
    if (c == 1) { out[0] = out[1] = out[2] = num(v, 0); return true; }
    if (c != 3) return false;
    for (int a = 0; a < 3; ++a) out[a] = num(v, a);
    return true;
}

// A box: a tuple (lo, hi) of three numbers each, or float[2][3]
static bool box_of(const vial_val_t* v, double lo[3], double hi[3]) {
    if (!v) return false;
    if (v->dtype == VIAL_DT_TUPLE) return v->num_elems == 2 && three(&v->elems[0], lo) && three(&v->elems[1], hi);
    if (!is_num(v) || vial_val_count(v) != 6) return false;
    for (int a = 0; a < 3; ++a) { lo[a] = num(v, a); hi[a] = num(v, 3 + a); }
    return true;
}

// ### SUMS ###

// An integer added from any thread: two's complement, so a negative share (a charge) adds as it should
static inline void add_share(int64_t* p, int64_t v) {
    md_atomic_fetch_add64((volatile uint64_t*)p, (uint64_t)v);
}

// w spread over the eight voxels around grid position g (voxel centres at whole numbers), by distance; corners outside
// the grid are dropped, or wrapped when it is periodic
static inline void deposit(int64_t* sum, const uint32_t dims[3], double gx, double gy, double gz, float w, bool periodic) {
    const double g[3] = {gx, gy, gz};
    int64_t i0[3];
    double f[3];
    for (int a = 0; a < 3; ++a) {
        const double fl = floor(g[a]);
        i0[a] = (int64_t)fl;
        f[a] = g[a] - fl;
    }
    for (int dx = 0; dx < 2; ++dx) {
        int64_t ix = i0[0] + dx;
        if (periodic) ix = ((ix % dims[0]) + dims[0]) % dims[0];
        else if (ix < 0 || ix >= dims[0]) continue;
        const double wx = dx ? f[0] : 1.0 - f[0];
        for (int dy = 0; dy < 2; ++dy) {
            int64_t iy = i0[1] + dy;
            if (periodic) iy = ((iy % dims[1]) + dims[1]) % dims[1];
            else if (iy < 0 || iy >= dims[1]) continue;
            const double wy = dy ? f[1] : 1.0 - f[1];
            for (int dz = 0; dz < 2; ++dz) {
                int64_t iz = i0[2] + dz;
                if (periodic) iz = ((iz % dims[2]) + dims[2]) % dims[2];
                else if (iz < 0 || iz >= dims[2]) continue;
                const double wz = dz ? f[2] : 1.0 - f[2];
                const int64_t c = (int64_t)llround((double)w * wx * wy * wz * VIAL_VOLUME_ONE);
                if (c) add_share(&sum[((size_t)ix * dims[1] + (size_t)iy) * dims[2] + (size_t)iz], c);
            }
        }
    }
}

static double voxel_volume(const vial_volume_t* vol) {
    const double* a = vol->basis[0];
    const double* b = vol->basis[1];
    const double* c = vol->basis[2];
    const double d = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
    return fabs(d);
}

static bool has_cell(const md_unitcell_t* cell) {
    return (md_unitcell_flags(cell) & (MD_UNITCELL_ORTHO | MD_UNITCELL_TRICLINIC)) != 0;
}

static double cell_volume(const md_unitcell_t* cell) {
    return has_cell(cell) ? fabs(cell->x * cell->y * cell->z) : 0.0;
}

// ### SETTING UP ###

// The positions and masses of particles idx[n] of xyz, as v1 weighs a structure
static void gather(vec4_t* out, const float* xyz, size_t num_atoms, const uint32_t* idx, uint32_t n, const md_system_t* sys) {
    for (uint32_t i = 0; i < n; ++i) {
        const uint32_t p = idx[i];
        const float* q = p < num_atoms ? xyz + 3 * (size_t)p : NULL;
        const float m = sys && p < sys->atom.count ? md_atom_mass(&sys->atom, p) : 1.0f;
        out[i] = q ? vec4_set(q[0], q[1], q[2], m > 0.0f ? m : 1.0f) : vec4_set(0, 0, 0, 0);
    }
}

bool vial_volume_begin(vial_kctx_t* k, vial_volume_t* vol, const char* name, const vial_val_t* const* args, const str_t* names,
                       uint32_t n, const vial_frame_t* first, bool add, md_allocator_i* alloc) {
    memset(vol, 0, sizeof(*vol));
    vol->alloc = alloc;
    const vargs_t a = {args, names, n};
    const bool sdf = !strcmp(name, "sdf");
    const vial_val_t* samples = sdf ? varg(&a, "samples", 3) : varg(&a, "samples", 1);
    const vial_val_t* sigma = sdf ? varg(&a, "sigma", 6) : varg(&a, "sigma", 4);
    const double want = is_num(samples) ? num(samples, 0) : VIAL_VOLUME_SAMPLES;
    if (!(want > 0.0)) return vial_rt_fail(k, "R0103", "samples= is a number of samples per Å, more than 0");
    vol->sigma = is_num(sigma) ? num(sigma, 0) : 0.0;
    if (!(vol->sigma >= 0.0)) return vial_rt_fail(k, "R0103", "sigma= is a length, 0 or more");
    double len[3];
    if (sdf) {
        vol->kind = VIAL_VOLUME_SDF;
        const vial_val_t* st = varg(&a, "structures", 0);
        const vial_val_t* ext = varg(&a, "extent", 2);
        const vial_val_t* ref = varg(&a, "ref", 4);
        if (!st || st->dtype != VIAL_DT_SET || !st->sets || !st->sets->count) return vial_rt_fail(k, "R0103", "sdf takes a set, or an array of sets, of structures");
        if (!three(ext, vol->extent) || !(vol->extent[0] > 0.0 && vol->extent[1] > 0.0 && vol->extent[2] > 0.0)) {
            return vial_rt_fail(k, "R0103", "extent= is a length more than 0, or three");
        }
        vol->structures = st->sets;
        vol->n = vial_sets_size(st->sets, 0);
        if (!vol->n) return vial_rt_fail(k, "R0103", "the structures of sdf are empty");
        for (uint32_t s = 1; s < st->sets->count; ++s) {
            if (vial_sets_size(st->sets, s) != vol->n) return vial_rt_fail(k, "R0103", "the structures of sdf have different sizes");
        }
        vol->ref_xyzw = md_alloc(alloc, vol->n * sizeof(vec4_t));
        const md_system_t* sys = k->in->msys;
        if (!is_none(ref)) {
            vial_val_t r = *ref;
            if (r.dtype == VIAL_DT_FLOAT64) r = vial_val_as(ref, VIAL_DT_FLOAT, k->alloc);
            if (r.dtype != VIAL_DT_FLOAT || vial_val_count(&r) != 3 * (size_t)vol->n || !r.data) {
                return vial_rt_fail(k, "R0103", "ref= has the positions of one structure: float[%u][3]", vol->n);
            }
            const uint32_t* idx = vial_sets_items(st->sets, 0);
            for (uint32_t i = 0; i < vol->n; ++i) {
                const float* q = (const float*)r.data + 3 * (size_t)i;
                const float m = sys && idx[i] < sys->atom.count ? md_atom_mass(&sys->atom, idx[i]) : 1.0f;
                vol->ref_xyzw[i] = vec4_set(q[0], q[1], q[2], m > 0.0f ? m : 1.0f);
            }
        } else {
            if (!first || !first->xyz) return vial_rt_fail(k, "R0104", "sdf fits onto frame 0, which the run does not have");
            gather(vol->ref_xyzw, first->xyz, first->num_atoms, vial_sets_items(st->sets, 0), vol->n, sys);
        }
        // Mutually consistent images about the reference's own centre, as v1 places it
        const md_unitcell_t cell = first ? first->cell : md_unitcell_none();
        md_util_deperiodize_self_vec4(vol->ref_xyzw, vol->n, &cell, &vol->ref_com);
        const mat3_eigen_t eig = mat3_eigen(mat3_covariance_matrix_vec4(vol->ref_xyzw, 0, vol->n, vol->ref_com));
        vol->axes = mat3_transpose(eig.vectors);
        for (int q = 0; q < 3; ++q) {
            len[q] = 2.0 * vol->extent[q];
            vol->origin[q] = -vol->extent[q];
        }
    } else {
        const vial_val_t* box = varg(&a, "box", 2);
        if (!is_none(box)) {
            vol->kind = VIAL_VOLUME_BOX;
            double hi[3];
            if (!box_of(box, vol->lo, hi)) return vial_rt_fail(k, "R0103", "box= is (lo, hi), each three lengths");
            for (int q = 0; q < 3; ++q) {
                len[q] = hi[q] - vol->lo[q];
                if (!(len[q] > 0.0)) return vial_rt_fail(k, "R0103", "box= has hi above lo along every axis");
                vol->len[q] = len[q];
                vol->origin[q] = vol->lo[q];
            }
        } else {
            vol->kind = VIAL_VOLUME_CELL;
            if (!first || !has_cell(&first->cell)) return vial_rt_fail(k, "R0103", "density over the unit cell needs a cell; give box=");
            float A[3][3];
            md_unitcell_A_extract_float(A, &first->cell);
            for (int q = 0; q < 3; ++q) len[q] = sqrt((double)A[q][0] * A[q][0] + (double)A[q][1] * A[q][1] + (double)A[q][2] * A[q][2]);
        }
    }
    vol->samples = vial_volume_grid(len, want, vol->dims, NULL);
    if (vol->kind == VIAL_VOLUME_CELL) {
        float A[3][3];
        md_unitcell_A_extract_float(A, &first->cell);
        for (int q = 0; q < 3; ++q) for (int c = 0; c < 3; ++c) vol->basis[q][c] = A[q][c] / vol->dims[q];
    } else {
        for (int q = 0; q < 3; ++q) vol->basis[q][q] = len[q] / vol->dims[q];
    }
    vol->vref = voxel_volume(vol);
    if (!add) return true;
    const size_t nv = (size_t)vol->dims[0] * vol->dims[1] * vol->dims[2];
    if (!vial_k_room(k, nv)) return false;
    vol->sum = md_alloc(md_get_heap_allocator(), nv * sizeof(int64_t));
    if (!vol->sum) return vial_rt_fail(k, "R0105", "a volume of %zu voxels does not fit in memory", nv);
    memset(vol->sum, 0, nv * sizeof(int64_t));
    vol->have_lock = md_mutex_init(&vol->lock);
    return true;
}

bool vial_volume_like(vial_volume_t* out, const vial_volume_t* vol, md_allocator_i* alloc) {
    *out = *vol;
    out->alloc = alloc;
    const size_t nv = (size_t)vol->dims[0] * vol->dims[1] * vol->dims[2];
    out->sum = md_alloc(md_get_heap_allocator(), nv * sizeof(int64_t));
    if (!out->sum) return false;
    memset(out->sum, 0, nv * sizeof(int64_t));
    if (vol->ref_xyzw) {
        out->ref_xyzw = md_alloc(alloc, vol->n * sizeof(vec4_t));
        memcpy(out->ref_xyzw, vol->ref_xyzw, vol->n * sizeof(vec4_t));
    }
    out->added = out->frames = out->bulk_frames = 0;
    out->bulk = 0.0;
    memset(out->cell_sum, 0, sizeof(out->cell_sum));
    out->have_lock = md_mutex_init(&out->lock);
    return true;
}

void vial_volume_free(vial_volume_t* vol) {
    if (!vol) return;
    if (vol->sum) md_free(md_get_heap_allocator(), vol->sum, (size_t)vol->dims[0] * vol->dims[1] * vol->dims[2] * sizeof(int64_t));
    if (vol->have_lock) md_mutex_destroy(&vol->lock);
    memset(vol, 0, sizeof(*vol));
}

double vial_volume_bulk(const vial_volume_t* vol) {
    return vol->bulk_frames ? vol->bulk / (double)vol->bulk_frames : NAN;
}

// ### ADDING A FRAME ###

// Structure s fitted at positions xyz: the rotation onto the reference and its centre, in the image it arrived in.
// False when it is too spread over the cell to place.
static bool fit_structure(const vial_volume_t* vol, uint32_t s, const float* xyz, size_t num_atoms, const md_unitcell_t* cell,
                          const md_system_t* sys, vec4_t* scratch, mat3_t* R, vec3_t* com) {
    gather(scratch, xyz, num_atoms, vial_sets_items(vol->structures, s), vol->n, sys);
    const float residual = md_util_optimal_rotation_pbc_vec4_iter(R, com, NULL, vol->ref_xyzw, vol->ref_com, scratch, vol->n, cell, 8, 1.0e-6f);
    // The residual is the margin to an ambiguous choice of images: near half a cell no choice per point represents
    // this copy, so it is left out rather than smeared into the volume in a wrong orientation
    float min_ext = FLT_MAX;
    if (has_cell(cell)) {
        float ext[3];
        md_unitcell_diag_extract_float(ext, cell);
        for (int q = 0; q < 3; ++q) if (ext[q] > 0.0f) min_ext = MIN(min_ext, ext[q]);
    }
    return residual <= 0.25f * min_ext;
}

bool vial_volume_pose(const vial_volume_t* vol, uint32_t s, const float* xyz, size_t num_atoms, const md_unitcell_t* cell, mat4_t* out) {
    if (!vol || vol->kind != VIAL_VOLUME_SDF || s >= vol->structures->count) return false;
    md_allocator_i* heap = md_get_heap_allocator();
    vec4_t* scratch = md_alloc(heap, vol->n * sizeof(vec4_t));
    mat3_t R;
    vec3_t com;
    const bool ok = fit_structure(vol, s, xyz, num_atoms, cell, NULL, scratch, &R, &com);
    md_free(heap, scratch, vol->n * sizeof(vec4_t));
    if (!ok) return false;
    const mat4_t RT = mat4_mul(mat4_from_mat3(R), mat4_translate(-com.x, -com.y, -com.z));
    *out = mat4_mul(mat4_from_mat3(vol->axes), RT);
    return true;
}

typedef struct sdf_cb_t {
    const vial_volume_t*    vol;
    mat3_t                  Q;          // the run's frame to the volume's
    vec3_t                  centre;     // where the points come in around
    const uint32_t*         own;        // the structure's particles, ascending, which are not its target
    uint32_t                num_own;
    const uint32_t*         ids;        // each point's particle, or NULL when not known
    const float*            w;          // each point's weight divided by a voxel's volume
    float                   w1;         // the weight of every point without w
    double                  inv_step[3];
} sdf_cb_t;

static bool owns(const uint32_t* own, uint32_t n, uint32_t p) {
    uint32_t lo = 0, hi = n;
    while (lo < hi) {
        const uint32_t mid = (lo + hi) / 2;
        if (own[mid] < p) lo = mid + 1;
        else hi = mid;
    }
    return lo < n && own[lo] == p;
}

static void sdf_cb(const uint32_t* idx, const float* x, const float* y, const float* z, size_t num_points, void* user) {
    const sdf_cb_t* c = user;
    const vial_volume_t* vol = c->vol;
    for (size_t i = 0; i < num_points; ++i) {
        const uint32_t pi = idx[i];
        if (c->ids && owns(c->own, c->num_own, c->ids[pi])) continue;
        const float w = c->w ? c->w[pi] : c->w1;
        if (w == 0.0f) continue;
        const vec3_t d = {{x[i] - c->centre.x, y[i] - c->centre.y, z[i] - c->centre.z}};
        const vec3_t l = mat3_mul_vec3(c->Q, d);
        // Grid positions, voxel centres at whole numbers
        const double gx = ((double)l.x - vol->origin[0]) * c->inv_step[0] - 0.5;
        const double gy = ((double)l.y - vol->origin[1]) * c->inv_step[1] - 0.5;
        const double gz = ((double)l.z - vol->origin[2]) * c->inv_step[2] - 0.5;
        if (gx < -1.0 || gy < -1.0 || gz < -1.0 || gx >= vol->dims[0] || gy >= vol->dims[1] || gz >= vol->dims[2]) continue;
        deposit(vol->sum, vol->dims, gx, gy, gz, w, false);
    }
}

bool vial_volume_add(vial_kctx_t* k, vial_volume_t* vol, const vial_val_t* const* args, const str_t* names, uint32_t n) {
    const vial_frame_t* fr = k->in->frame;
    if (!fr || !fr->xyz) return vial_rt_fail(k, "R0103", "a volume adds frames");
    const vargs_t a = {args, names, n};
    const bool sdf = vol->kind == VIAL_VOLUME_SDF;
    const vial_val_t* tv = sdf ? varg(&a, "target", 1) : varg(&a, "points", 0);
    const vial_val_t* wv = sdf ? varg(&a, "weights", 5) : varg(&a, "weights", 3);
    if (!tv) return vial_rt_fail(k, "R0103", "a missing argument");
    vial_val_t p;
    if (!vial_k_positions(k, tv, &p)) return false;
    if (p.dtype == VIAL_DT_FLOAT64) p = vial_val_as(&p, VIAL_DT_FLOAT, k->alloc);
    if (p.dtype != VIAL_DT_FLOAT || !p.ndim || p.dims[p.ndim - 1] != 3) return vial_rt_fail(k, "R0103", "a volume takes points, float[n][3], or a set");
    const size_t np = vial_val_count(&p) / 3;
    // Which particle each point is: a set's, or what the points are positions of
    const uint32_t* ids = NULL;
    if (tv->dtype == VIAL_DT_SET && tv->sets && tv->sets->count == 1 && vial_sets_size(tv->sets, 0) == np) ids = vial_sets_items(tv->sets, 0);
    else if (p.origin && p.origin->count == 1 && vial_sets_size(p.origin, 0) == np) ids = vial_sets_items(p.origin, 0);
    // Weights, over a cell that changes scaled to frame 0's voxels; a missing point weighs nothing
    const double vox = vol->kind == VIAL_VOLUME_CELL ? cell_volume(&fr->cell) / ((double)vol->dims[0] * vol->dims[1] * vol->dims[2]) : vol->vref;
    if (vol->kind == VIAL_VOLUME_CELL && has_cell(&fr->cell) && !(vox > 0.0)) return vial_rt_fail(k, "R0103", "a volume over a cell of no volume");
    const double scale = vox > 0.0 ? vol->vref / vox : 1.0;
    float* w = md_alloc(k->alloc, (np ? np : 1) * sizeof(float));
    double total = 0.0;
    if (!is_none(wv)) {
        if (!is_num(wv) || vial_val_count(wv) != np) return vial_rt_fail(k, "R0103", "weights= has one number per point: %zu, not %zu", np, wv ? vial_val_count(wv) : 0);
        for (size_t i = 0; i < np; ++i) {
            const double x = vial_val_is_missing(wv, i) ? 0.0 : num(wv, i);
            w[i] = (float)(x * scale);
            total += x;
        }
    } else {
        for (size_t i = 0; i < np; ++i) w[i] = (float)scale;
        total = (double)np;
    }
    for (size_t i = 0; i < np; ++i) {
        const float* q = (const float*)p.data + 3 * i;
        if (vial_val_is_missing(&p, 3 * i) || !isfinite(q[0]) || !isfinite(q[1]) || !isfinite(q[2])) w[i] = 0.0f;
    }
    const md_unitcell_t cell = fr->cell;
    uint64_t added = 0;
    if (sdf) {
        const uint32_t ns = vol->structures->count;
        vec4_t* scratch = md_alloc(k->alloc, vol->n * sizeof(vec4_t));
        double reach = 0.0;
        for (int q = 0; q < 3; ++q) reach = MAX(reach, vol->extent[q] + sqrt(vol->basis[q][0] * vol->basis[q][0] + vol->basis[q][1] * vol->basis[q][1] + vol->basis[q][2] * vol->basis[q][2]));
        md_spatial_acc_t acc = {.alloc = k->alloc};
        md_coord_stream_t pts = md_coord_stream_from_aos((const float*)p.data, 3 * sizeof(float), NULL, np);
        // The cell extent of the last frame this interpreter built for, so the spread of the points is estimated once
        vial_interp_t* in = k->in;
        const size_t opi = k->op ? (size_t)(k->op - in->ir->ops) : SIZE_MAX;
        double* ext = in->cell_ext && opi < md_array_size(in->cell_ext) ? &in->cell_ext[opi] : NULL;
        if (np) {
            md_spatial_acc_init(&acc, &(md_spatial_acc_desc_t){.coords = &pts, .cutoff = reach, .unitcell = &cell, .cell_ext = ext ? *ext : 0.0});
            if (ext) *ext = md_spatial_acc_cell_ext(&acc);
        }
        for (uint32_t s = 0; s < ns; ++s) {
            mat3_t R;
            vec3_t com;
            if (!fit_structure(vol, s, fr->xyz, fr->num_atoms, &cell, k->in->msys, scratch, &R, &com)) continue;
            added += 1;
            if (!np) continue;
            // The points come in around the query's centre, folded into the cell: a lattice vector off com when the
            // structure sits outside it
            const double cen[3] = {com.x, com.y, com.z};
            double qc[3];
            md_spatial_acc_aabb_query_center(qc, &acc, cen);
            sdf_cb_t c = {.vol = vol, .Q = mat3_mul(vol->axes, R), .centre = {{(float)qc[0], (float)qc[1], (float)qc[2]}},
                          .own = vial_sets_items(vol->structures, s), .num_own = vol->n, .ids = ids, .w = w};
            for (int q = 0; q < 3; ++q) c.inv_step[q] = vol->dims[q] / (2.0 * vol->extent[q]);
            // The box about it along each axis of the run's frame: what encloses the turned box, a voxel more for the
            // eight voxels a point spreads over
            double rad[3];
            for (int q = 0; q < 3; ++q) {
                rad[q] = 0.0;
                for (int j = 0; j < 3; ++j) rad[q] += fabs(c.Q.elem[q][j]) * (vol->extent[j] + vol->basis[j][j]);
            }
            md_spatial_acc_for_each_point_in_aabb(&acc, cen, rad, sdf_cb, &c);
        }
        if (np) md_spatial_acc_free(&acc);
    } else if (vol->kind == VIAL_VOLUME_BOX) {
        const bool periodic = has_cell(&cell);
        double inv[3], c[3];
        for (int q = 0; q < 3; ++q) { inv[q] = vol->dims[q] / vol->len[q]; c[q] = vol->lo[q] + 0.5 * vol->len[q]; }
        for (size_t i = 0; i < np; ++i) {
            if (w[i] == 0.0f) continue;
            const float* q = (const float*)p.data + 3 * i;
            // The image nearest the box's centre
            vec3_t d = {{(float)(q[0] - c[0]), (float)(q[1] - c[1]), (float)(q[2] - c[2])}};
            if (periodic) md_util_min_image_vec3(&d, 1, &cell);
            const double gx = (c[0] + d.x - vol->lo[0]) * inv[0] - 0.5;
            const double gy = (c[1] + d.y - vol->lo[1]) * inv[1] - 0.5;
            const double gz = (c[2] + d.z - vol->lo[2]) * inv[2] - 0.5;
            if (gx < -1.0 || gy < -1.0 || gz < -1.0 || gx >= vol->dims[0] || gy >= vol->dims[1] || gz >= vol->dims[2]) continue;
            deposit(vol->sum, vol->dims, gx, gy, gz, w[i], false);
        }
        added = 1;
    } else {
        if (!has_cell(&cell)) {
            vial_rt_note(k, "R0101", "frame %lld has no unit cell; a density over the cell leaves it out", (long long)fr->index);
            return true;
        }
        float I[3][3];
        md_unitcell_I_extract_float(I, &cell);
        for (size_t i = 0; i < np; ++i) {
            if (w[i] == 0.0f) continue;
            const float* q = (const float*)p.data + 3 * i;
            double g[3];
            for (int r = 0; r < 3; ++r) {
                // Fractional, I[col][row] as A is
                double s = (double)I[0][r] * q[0] + (double)I[1][r] * q[1] + (double)I[2][r] * q[2];
                s -= floor(s);
                g[r] = s * vol->dims[r] - 0.5;
            }
            deposit(vol->sum, vol->dims, g[0], g[1], g[2], w[i], true);
        }
        added = 1;
    }
    const double cv = cell_volume(&cell);
    if (vol->have_lock) md_mutex_lock(&vol->lock);
    vol->added += added;
    vol->frames += 1;
    if (cv > 0.0) { vol->bulk += total / cv; vol->bulk_frames += 1; }
    if (vol->kind == VIAL_VOLUME_CELL && cv > 0.0) {
        // The mean cell, which the published grid is laid out in
        float A[3][3];
        md_unitcell_A_extract_float(A, &cell);
        for (int q = 0; q < 3; ++q) for (int c = 0; c < 3; ++c) vol->cell_sum[q][c] += A[q][c];
    }
    if (vol->have_lock) md_mutex_unlock(&vol->lock);
    return true;
}

// ### THE END ###

// Each line along axis a convolved with a Gaussian of sigma s voxels, normalized over its whole width, in place: zero
// past the ends, or around them when periodic
static void smooth_axis(float* v, const uint32_t dims[3], int a, double s, bool periodic) {
    const int r = (int)ceil(3.0 * s);
    if (r < 1) return;
    float* kern = md_alloc(md_get_heap_allocator(), (2 * (size_t)r + 1) * sizeof(float));
    double tot = 0.0;
    for (int i = -r; i <= r; ++i) tot += exp(-0.5 * (double)i * i / (s * s));
    for (int i = -r; i <= r; ++i) kern[i + r] = (float)(exp(-0.5 * (double)i * i / (s * s)) / tot);
    const uint32_t n = dims[a];
    float* line = md_alloc(md_get_heap_allocator(), (size_t)n * sizeof(float));
    const size_t stride[3] = {(size_t)dims[1] * dims[2], dims[2], 1};
    const int b = (a + 1) % 3, c = (a + 2) % 3;
    for (uint32_t i = 0; i < dims[b]; ++i) {
        for (uint32_t j = 0; j < dims[c]; ++j) {
            float* base = v + i * stride[b] + j * stride[c];
            for (uint32_t x = 0; x < n; ++x) line[x] = base[x * stride[a]];
            for (uint32_t x = 0; x < n; ++x) {
                double acc = 0.0;
                for (int t = -r; t <= r; ++t) {
                    int64_t y = (int64_t)x + t;
                    if (periodic) y = ((y % n) + n) % n;
                    else if (y < 0 || y >= n) continue;
                    acc += (double)kern[t + r] * line[y];
                }
                base[x * stride[a]] = (float)acc;
            }
        }
    }
    md_free(md_get_heap_allocator(), line, (size_t)n * sizeof(float));
    md_free(md_get_heap_allocator(), kern, (2 * (size_t)r + 1) * sizeof(float));
}

bool vial_volume_end(vial_kctx_t* k, vial_volume_t* vol, vial_val_t* out) {
    (void)k;
    vial_volume_t* v = vol;
    const size_t nv = (size_t)v->dims[0] * v->dims[1] * v->dims[2];
    // Per Å³ and per what was added, as floats in place of the integers: float i is written over the first half of
    // integer i/2, which was read before
    const double scale = v->added && v->vref > 0.0 ? 1.0 / (VIAL_VOLUME_ONE * v->vref * (double)v->added) : 0.0;
    uint8_t* bytes = (uint8_t*)v->sum;
    for (size_t i = 0; i < nv; ++i) {
        int64_t s;
        memcpy(&s, bytes + 8 * i, 8);
        const float x = (float)((double)s * scale);
        memcpy(bytes + 4 * i, &x, 4);
    }
    float* f = (float*)v->sum;
    if (v->kind == VIAL_VOLUME_CELL && v->bulk_frames) {
        for (int q = 0; q < 3; ++q) for (int c = 0; c < 3; ++c) v->basis[q][c] = v->cell_sum[q][c] / (double)v->bulk_frames / v->dims[q];
    }
    if (v->sigma > 0.0) {
        for (int a = 0; a < 3; ++a) {
            const double step = sqrt(v->basis[a][0] * v->basis[a][0] + v->basis[a][1] * v->basis[a][1] + v->basis[a][2] * v->basis[a][2]);
            if (step > 0.0) smooth_axis(f, v->dims, a, v->sigma / step, v->kind == VIAL_VOLUME_CELL);
        }
    }
    *out = (vial_val_t){.dtype = VIAL_DT_FLOAT, .ndim = 3, .dims = {v->dims[0], v->dims[1], v->dims[2]}, .data = f};
    return true;
}

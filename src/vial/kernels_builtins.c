// The interpreter's built-ins: geometry, properties, reductions, math and the array functions. Each takes its
// arguments by parameter, as the check bound them. Built with floating-point contraction off, like kernels.c.

// No contraction: where a multiply and an add fuse is the IR's decision (VIAL_IR_FMA). GCC and Clang get
// -ffp-contract=off for this file from CMake; MSVC's /fp:contract is taken back here, which it allows per file.
#if defined(_MSC_VER) && !defined(__clang__)
#pragma fp_contract(off)
#endif

#include "interp.h"
#include "infer_internal.h"
#include "ast.h"

#include <md_system.h>
#include <md_util.h>
#include <core/md_intrinsics.h>
#include <core/md_simd.h>
#include <core/md_spatial_acc.h>
#include <core/md_coord_stream.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_vec_math.h>

#include <float.h>
#include <math.h>
#include <string.h>

#define DEG_PER_RAD 57.29577951308232
#define RAD_PER_DEG 0.017453292519943295

bool vial_k_index_set(vial_kctx_t* k, const vial_val_t* a, const vial_index_item_t* it, vial_val_t* out);

typedef struct bargs_t {
    const vial_val_t* const* v;
    const str_t*      names;
    uint32_t          n;
} bargs_t;

// An argument by parameter name when the check bound them, or else the pos-th positional one; NULL when not given
static const vial_val_t* arg(const bargs_t* a, const char* name, uint32_t pos) {
    bool bound = false;
    for (uint32_t i = 0; i < a->n; ++i) {
        if (!a->names[i].len) continue;
        if (a->names[i].ptr[0] == name[0] && str_eq(a->names[i], str_from_cstr(name))) return a->v[i];
        if (i == 0) bound = true;   // a bound call names its first parameter
    }
    if (bound) return NULL;
    uint32_t p = 0;
    for (uint32_t i = 0; i < a->n; ++i) {
        if (a->names[i].len) continue;
        if (p++ == pos) return a->v[i];
    }
    return NULL;
}

// The type of the i-th argument as the IR has it
static const vial_type_t* arg_type(vial_kctx_t* k, const bargs_t* a, const vial_val_t* v) {
    for (uint32_t i = 0; i < a->n; ++i) {
        if (a->v[i] == v) return &k->in->ir->values[k->in->ir->args[k->op->arg_beg + i]].type;
    }
    return NULL;
}

static double num_at(const vial_val_t* v, size_t i) {
    switch (v->dtype) {
    case VIAL_DT_BOOL:    return ((const uint8_t*)v->data)[i];
    case VIAL_DT_INT:     return ((const int32_t*)v->data)[i];
    case VIAL_DT_FLOAT:   return ((const float*)v->data)[i];
    case VIAL_DT_FLOAT64: return ((const double*)v->data)[i];
    default:              return 0.0;
    }
}

static bool is_num(uint8_t dt) { return dt == VIAL_DT_BOOL || dt == VIAL_DT_INT || dt == VIAL_DT_FLOAT || dt == VIAL_DT_FLOAT64; }


static bool str_is(const vial_val_t* v, const char* s) {
    return v && v->dtype == VIAL_DT_STR && v->ndim == 0 && str_eq(*(const str_t*)v->data, str_from_cstr(s));
}

// ### POSITIONS ###

// Each particle's mass, looked up once for the interpreter's life, as v1 did: a property of the system, not of a frame
static const float* k_masses(vial_kctx_t* k) {
    vial_interp_t* in = k->in;
    if (in->masses || !in->msys) return in->masses;
    const size_t n = in->msys->atom.count;
    float* m = md_alloc(in->alloc, (n ? n : 1) * sizeof(float));
    for (size_t i = 0; i < n; ++i) m[i] = md_atom_mass(&in->msys->atom, i);
    in->masses = m;
    return m;
}

static inline float k_mass(const float* masses, const md_system_t* sys, size_t i) {
    return i < sys->atom.count ? masses[i] : 0.0f;
}

static const float* frame_xyz(vial_kctx_t* k) {
    if (!k->in->frame || !k->in->frame->xyz) { vial_rt_fail(k, "R0103", "positions need a frame"); return NULL; }
    return k->in->frame->xyz;
}

// The positions of particles idx[n] into out[3n]: one copy where they are a run, as slices and residues are (a set's
// particles are sorted and distinct, so they are a run when the last is n - 1 past the first). False past the frame.
static bool gather_xyz(float* out, const float* xyz, size_t na, const uint32_t* idx, uint32_t n) {
    if (!n) return true;
    if (idx[n - 1] >= na || idx[0] >= na) return false;
    if (idx[n - 1] - idx[0] == n - 1) {
        memcpy(out, xyz + 3 * (size_t)idx[0], 3 * sizeof(float) * n);
        return true;
    }
    for (uint32_t i = 0; i < n; ++i) {
        if (idx[i] >= na) return false;
        const float* p = xyz + 3 * (size_t)idx[i];
        out[3 * i + 0] = p[0];
        out[3 * i + 1] = p[1];
        out[3 * i + 2] = p[2];
    }
    return true;
}

bool vial_k_positions(vial_kctx_t* k, const vial_val_t* s, vial_val_t* out) {
    if (!s) return vial_rt_fail(k, "R0103", "a missing argument");
    if (s->dtype != VIAL_DT_SET) { *out = *s; return true; }   // coordinates are themselves
    const float* xyz = frame_xyz(k);
    if (!xyz) return false;
    const size_t na = k->in->frame->num_atoms;
    const vial_sets_t* set = s->sets;
    if (s->ndim == 0) {
        const uint32_t n = vial_sets_size(set, 0);
        vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, 2, (uint32_t[]){n, 3}, &o)) return false;
        const uint32_t* idx = vial_sets_items(set, 0);
        if (!gather_xyz(o.data, xyz, na, idx, n)) return vial_rt_fail(k, "R0103", "a particle is beyond the frame's %zu", na);
        o.origin = set;
        *out = o;
        return true;
    }
    // An array of sets of one particle each: one position per set
    const size_t cnt = vial_val_count(s);
    uint32_t dims[VIAL_MAX_DIMS + 1];
    memcpy(dims, s->dims, s->ndim * sizeof(uint32_t));
    dims[s->ndim] = 3;
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, s->ndim + 1, dims, &o)) return false;
    for (size_t i = 0; i < cnt; ++i) {
        const uint32_t sz = vial_sets_size(set, (uint32_t)i);
        if (sz == 1) {
            const uint32_t a = vial_sets_items(set, (uint32_t)i)[0];
            if (a >= na) return vial_rt_fail(k, "R0103", "particle %u is beyond the frame's %zu", a, na);
            memcpy((float*)o.data + 3 * i, xyz + 3 * (size_t)a, 3 * sizeof(float));
        } else if (sz == 0) {
            for (int c = 0; c < 3; ++c) vial_val_set_missing(&o, 3 * i + c, k->alloc);
        } else {
            return vial_rt_fail(k, "R0103", "positions of sets of several particles do not form one array");
        }
    }
    *out = o;
    return true;
}

// A spatial search of this op's: built with the cell extent its last one was, so the spread of the points is
// estimated once, in the first frame, not in every frame (md_spatial_acc_desc_t.cell_ext)
static void k_spatial_init(vial_kctx_t* k, md_spatial_acc_t* acc, md_coord_stream_t* coords, double cutoff, const md_unitcell_t* cell) {
    vial_interp_t* in = k->in;
    double* ext = NULL;
    if (k->op && in->cell_ext) {
        const size_t opi = (size_t)(k->op - in->ir->ops);
        if (opi < md_array_size(in->cell_ext)) ext = &in->cell_ext[opi];
    }
    md_spatial_acc_init(acc, &(md_spatial_acc_desc_t){.coords = coords, .cutoff = cutoff, .unitcell = cell, .cell_ext = ext ? *ext : 0.0});
    if (ext) *ext = md_spatial_acc_cell_ext(acc);
}

static md_unitcell_t frame_cell(vial_kctx_t* k) {
    return k->in->frame ? k->in->frame->cell : (md_unitcell_t){0};
}

// ### GEOMETRY ###

// Points from a value: positions of a set, or coordinates; the last dimension is 3
static bool points(vial_kctx_t* k, const vial_val_t* v, vial_val_t* out) {
    if (!vial_k_positions(k, v, out)) return false;
    if (out->dtype != VIAL_DT_FLOAT) *out = vial_val_as(out, VIAL_DT_FLOAT, k->alloc);
    if (!out->ndim || out->dims[out->ndim - 1] != 3) return vial_rt_fail(k, "R0103", "points have a last dimension of 3");
    return true;
}

// A point parameter's argument: a set of one particle is its position, float[3]; an array of them, one per set
static bool point(vial_kctx_t* k, const vial_val_t* v, vial_val_t* out) {
    if (v->dtype == VIAL_DT_SET && v->ndim == 0) return vial_k_kind(k, VX_KIND_POINT, v, out);
    return points(k, v, out);
}

// b - a under the minimum image, broadcast over leading dimensions
static bool k_displacement(vial_kctx_t* k, const vial_val_t* a, const vial_val_t* b, vial_val_t* out) {
    if (!a || !b) return vial_rt_fail(k, "R0103", "a missing argument");
    vial_val_t pa, pb;
    if (!point(k, a, &pa) || !point(k, b, &pb)) return false;
    // Broadcast the leading dimensions
    vial_val_t la = pa, lb = pb;
    la.ndim -= 1; lb.ndim -= 1;
    uint8_t nd = la.ndim > lb.ndim ? la.ndim : lb.ndim;
    uint32_t dims[VIAL_MAX_DIMS];
    for (uint8_t d = 0; d < nd; ++d) dims[d] = 1;
    for (int w = 0; w < 2; ++w) {
        const vial_val_t* v = w ? &lb : &la;
        for (uint8_t d = 0; d < v->ndim; ++d) {
            uint32_t* o = &dims[nd - v->ndim + d];
            if (*o == 1) *o = v->dims[d];
            else if (v->dims[d] != 1 && v->dims[d] != *o) return vial_rt_fail(k, "R0103", "extents %u and %u do not broadcast", *o, v->dims[d]);
        }
    }
    size_t n = 1;
    for (uint8_t d = 0; d < nd; ++d) n *= dims[d];
    uint32_t odims[VIAL_MAX_DIMS + 1];
    memcpy(odims, dims, nd * sizeof(uint32_t));
    odims[nd] = 3;
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, nd + 1, odims, &o)) return false;
    const md_unitcell_t cell = frame_cell(k);
    for (size_t i = 0; i < n; ++i) {
        size_t rem = i, ia = 0, ib = 0, sa = 1, sb = 1;
        for (int d = nd - 1; d >= 0; --d) {
            const uint32_t x = (uint32_t)(rem % dims[d]);
            rem /= dims[d];
            const int da = d - (nd - la.ndim), db = d - (nd - lb.ndim);
            if (da >= 0) { if (la.dims[da] != 1) ia += x * sa; sa *= la.dims[da]; }
            if (db >= 0) { if (lb.dims[db] != 1) ib += x * sb; sb *= lb.dims[db]; }
        }
        const float* A = (const float*)pa.data + 3 * ia;
        const float* B = (const float*)pb.data + 3 * ib;
        float* O = (float*)o.data + 3 * i;
        O[0] = B[0] - A[0]; O[1] = B[1] - A[1]; O[2] = B[2] - A[2];
        bool miss = false;
        for (int c = 0; c < 3; ++c) miss = miss || vial_val_is_missing(&pa, 3 * ia + c) || vial_val_is_missing(&pb, 3 * ib + c);
        if (miss) for (int c = 0; c < 3; ++c) vial_val_set_missing(&o, 3 * i + c, k->alloc);
    }
    // The minimum image over every difference in one call
    md_util_min_image_vec3((vec3_t*)o.data, n, &cell);
    if (pb.origin && nd == lb.ndim) o.origin = pb.origin;
    *out = o;
    return true;
}

// Memory aligned to 64 bytes. mdlib's packed kernels peel scalar iterations up to the alignment of their input, so
// where an array starts decides how its sum is grouped: given a 64-aligned array, they group the same way wherever
// it is, and a result does not depend on the address of its operands.
static void* k_alloc64(vial_kctx_t* k, size_t bytes) {
    uint8_t* p = md_alloc(k->alloc, bytes + 64);
    return (void*)(((uintptr_t)p + 63) & ~(uintptr_t)63);
}

// The centre of mass of each set or group of points; periodic-aware
static bool k_com(vial_kctx_t* k, const vial_val_t* s, const vial_val_t* weights, vial_val_t* out) {
    if (!s) return vial_rt_fail(k, "R0103", "a missing argument");
    const bool by_mass = !weights || str_is(weights, "mass");
    const md_unitcell_t cell = frame_cell(k);
    const md_system_t* sys = k->in->msys;
    if (s->dtype == VIAL_DT_SET) {
        const float* xyz = frame_xyz(k);
        if (!xyz) return false;
        const size_t cnt = s->ndim ? vial_val_count(s) : 1;
        uint32_t dims[VIAL_MAX_DIMS + 1];
        memcpy(dims, s->dims, s->ndim * sizeof(uint32_t));
        dims[s->ndim] = 3;
        vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, s->ndim + 1, dims, &o)) return false;
        for (size_t i = 0; i < cnt; ++i) {
            const uint32_t n = vial_sets_size(s->sets, (uint32_t)i);
            const int32_t* idx = (const int32_t*)vial_sets_items(s->sets, (uint32_t)i);
            if (!n) { for (int c = 0; c < 3; ++c) vial_val_set_missing(&o, 3 * i + c, k->alloc); continue; }
            // Gathered, so mdlib's packed kernel runs over contiguous positions and weights
            vec3_t* g = k_alloc64(k, n * sizeof(vec3_t));
            float* w = by_mass && sys ? md_alloc(k->alloc, n * sizeof(float)) : NULL;
            if (!gather_xyz((float*)g, xyz, k->in->frame->num_atoms, (const uint32_t*)idx, n)) {
                return vial_rt_fail(k, "R0103", "a particle is beyond the frame's %zu", k->in->frame->num_atoms);
            }
            if (w) {
                const float* ms = k_masses(k);
                for (uint32_t q = 0; q < n; ++q) w[q] = k_mass(ms, sys, (size_t)idx[q]);
            }
            const vec3_t c = md_util_com_compute(g, w, NULL, n, &cell);
            float* O = (float*)o.data + 3 * i;
            O[0] = c.x; O[1] = c.y; O[2] = c.z;
        }
        *out = o;
        return true;
    }
    vial_val_t p;
    if (!points(k, s, &p)) return false;
    if (p.ndim < 2) { *out = p; return true; }
    const uint32_t n = p.dims[p.ndim - 2];
    const size_t groups = vial_val_count(&p) / (3 * (size_t)(n ? n : 1));
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, p.ndim - 1, p.dims, &o)) return false;
    o.dims[p.ndim - 2] = 3;
    float* w = NULL;
    if (by_mass) {
        if (!p.origin || !sys || vial_sets_size(p.origin, 0) != n) return vial_rt_fail(k, "R0103", "com by mass of points whose particles are not known; write weights=None");
        w = md_alloc(k->alloc, n * sizeof(float));
        const float* ms = k_masses(k);
        for (uint32_t q = 0; q < n; ++q) w[q] = k_mass(ms, sys, vial_sets_items(p.origin, 0)[q]);
    }
    // Each group copied to one aligned buffer, for the same grouping of its sum (k_alloc64)
    vec3_t* buf = n ? k_alloc64(k, n * sizeof(vec3_t)) : NULL;
    for (size_t g = 0; g < groups; ++g) {
        float* O = (float*)o.data + 3 * g;
        if (!n) { for (int c = 0; c < 3; ++c) vial_val_set_missing(&o, 3 * g + c, k->alloc); continue; }
        memcpy(buf, (const float*)p.data + 3 * (size_t)n * g, n * sizeof(vec3_t));
        const vec3_t c = md_util_com_compute(buf, w, NULL, n, &cell);
        O[0] = c.x; O[1] = c.y; O[2] = c.z;
    }
    *out = o;
    return true;
}

// The angle at b between a and c, in degrees
static float angle3(const float* a, const float* b, const float* c, const md_unitcell_t* cell) {
    vec3_t u = {{a[0] - b[0], a[1] - b[1], a[2] - b[2]}};
    vec3_t v = {{c[0] - b[0], c[1] - b[1], c[2] - b[2]}};
    md_util_min_image_vec3(&u, 1, cell);
    md_util_min_image_vec3(&v, 1, cell);
    const double uu = (double)u.x * u.x + (double)u.y * u.y + (double)u.z * u.z;
    const double vv = (double)v.x * v.x + (double)v.y * v.y + (double)v.z * v.z;
    const double uv = (double)u.x * v.x + (double)u.y * v.y + (double)u.z * v.z;
    double cs = uv / sqrt(uu * vv);
    if (cs > 1.0) cs = 1.0;
    if (cs < -1.0) cs = -1.0;
    return (float)(acos(cs) * DEG_PER_RAD);
}

// The dihedral of a-b-c-d, in degrees, in (-180, 180]
static float dihedral4(const float* a, const float* b, const float* c, const float* d, const md_unitcell_t* cell) {
    vec3_t b1 = {{b[0] - a[0], b[1] - a[1], b[2] - a[2]}};
    vec3_t b2 = {{c[0] - b[0], c[1] - b[1], c[2] - b[2]}};
    vec3_t b3 = {{d[0] - c[0], d[1] - c[1], d[2] - c[2]}};
    md_util_min_image_vec3(&b1, 1, cell);
    md_util_min_image_vec3(&b2, 1, cell);
    md_util_min_image_vec3(&b3, 1, cell);
    // Blondel and Karplus: atan2(|b2| b1 . (b2 x b3), (b1 x b2) . (b2 x b3)), positive clockwise, as IUPAC
    const double n1[3] = {(double)b1.y * b2.z - (double)b1.z * b2.y, (double)b1.z * b2.x - (double)b1.x * b2.z, (double)b1.x * b2.y - (double)b1.y * b2.x};
    const double n2[3] = {(double)b2.y * b3.z - (double)b2.z * b3.y, (double)b2.z * b3.x - (double)b2.x * b3.z, (double)b2.x * b3.y - (double)b2.y * b3.x};
    const double lb2 = sqrt((double)b2.x * b2.x + (double)b2.y * b2.y + (double)b2.z * b2.z);
    const double y = lb2 * ((double)b1.x * n2[0] + (double)b1.y * n2[1] + (double)b1.z * n2[2]);
    const double x = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
    double r = atan2(y, x) * DEG_PER_RAD;
    if (r <= -180.0) r += 360.0;
    return (float)r;
}

// angle(a, b, c) or angle(m) with m's last two dimensions [3][3]; dihedral likewise with 4
static bool k_angle(vial_kctx_t* k, const bargs_t* ba, int npts, vial_val_t* out) {
    const char* pn[4] = {"a", "b", "c", "d"};
    vial_val_t p[4];
    const vial_val_t* one = arg(ba, "points", 0);
    const vial_val_t* first = one ? one : arg(ba, "a", 0);
    const vial_val_t* second = one ? NULL : arg(ba, "b", 1);
    const md_unitcell_t cell = frame_cell(k);
    if (first && !second) {
        // One array of npts points as rows
        vial_val_t m;
        if (!points(k, first, &m)) return false;
        if (m.ndim < 2 || m.dims[m.ndim - 2] != (uint32_t)npts) return vial_rt_fail(k, "R0103", "%d points as rows", npts);
        const size_t g = vial_val_count(&m) / (3 * (size_t)npts);
        vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, m.ndim - 2, m.dims, &o)) return false;
        for (size_t i = 0; i < g; ++i) {
            const float* r = (const float*)m.data + 3 * npts * i;
            ((float*)o.data)[i] = npts == 3 ? angle3(r, r + 3, r + 6, &cell) : dihedral4(r, r + 3, r + 6, r + 9, &cell);
            for (int q = 0; q < 3 * npts; ++q) if (vial_val_is_missing(&m, 3 * npts * i + q)) { vial_val_set_missing(&o, i, k->alloc); break; }
        }
        *out = o;
        return true;
    }
    for (int i = 0; i < npts; ++i) {
        const vial_val_t* v = arg(ba, pn[i], (uint32_t)i);
        if (!v) return vial_rt_fail(k, "R0103", "%s needs %d points", npts == 3 ? "angle" : "dihedral", npts);
        if (!point(k, v, &p[i])) return false;
    }
    // Leading dimensions broadcast
    uint8_t nd = 0;
    for (int i = 0; i < npts; ++i) if (p[i].ndim - 1 > nd) nd = p[i].ndim - 1;
    uint32_t dims[VIAL_MAX_DIMS];
    for (uint8_t d = 0; d < nd; ++d) dims[d] = 1;
    for (int i = 0; i < npts; ++i) {
        const uint8_t pn2 = p[i].ndim - 1;
        for (uint8_t d = 0; d < pn2; ++d) {
            uint32_t* o = &dims[nd - pn2 + d];
            if (*o == 1) *o = p[i].dims[d];
            else if (p[i].dims[d] != 1 && p[i].dims[d] != *o) return vial_rt_fail(k, "R0103", "extents do not broadcast");
        }
    }
    size_t n = 1;
    for (uint8_t d = 0; d < nd; ++d) n *= dims[d];
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, nd, dims, &o)) return false;
    for (size_t i = 0; i < n; ++i) {
        const float* q[4];
        bool miss = false;
        for (int j = 0; j < npts; ++j) {
            const uint8_t pn2 = p[j].ndim - 1;
            size_t rem = i, off = 0, s = 1;
            for (int d = nd - 1; d >= 0; --d) {
                const uint32_t x = (uint32_t)(rem % dims[d]);
                rem /= dims[d];
                const int dd = d - (nd - pn2);
                if (dd >= 0) { if (p[j].dims[dd] != 1) off += x * s; s *= p[j].dims[dd]; }
            }
            q[j] = (const float*)p[j].data + 3 * off;
            for (int c = 0; c < 3; ++c) miss = miss || vial_val_is_missing(&p[j], 3 * off + c);
        }
        ((float*)o.data)[i] = npts == 3 ? angle3(q[0], q[1], q[2], &cell) : dihedral4(q[0], q[1], q[2], q[3], &cell);
        if (miss) vial_val_set_missing(&o, i, k->alloc);
    }
    *out = o;
    return true;
}

// The closest or farthest pair of points of p and q under the minimum image: (distance, i, j). Distinct leaves out
// the pairs of a particle with itself, by the points' particles; or, mapped over a dimension that runs over the units
// q's points do, the pair of the call's own unit, skip.
static bool k_pair_x(vial_kctx_t* k, const vial_val_t* pv, const vial_val_t* qv, bool closest, bool distinct, int64_t skip, vial_val_t* out);
static bool k_pair(vial_kctx_t* k, const vial_val_t* pv, const vial_val_t* qv, bool closest, vial_val_t* out) {
    return k_pair_x(k, pv, qv, closest, false, -1, out);
}

// The matrix's entry for points a and b: the length of b - a under the minimum image, as distance() computes it
static inline float pair_entry(const float* a, const float* b, const md_unitcell_t* cell) {
    vec3_t v = {{b[0] - a[0], b[1] - a[1], b[2] - a[2]}};
    md_util_min_image_vec3(&v, 1, cell);
    return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

// Row i of the matrix, each entry as pair_entry gives it, into d[m]; row holds m
static void pair_row(const float* a, const float* Q, uint32_t m, const md_unitcell_t* cell, vec3_t* row, float* d) {
    for (uint32_t j = 0; j < m; ++j) {
        const float* b = Q + 3 * j;
        row[j] = (vec3_t){{b[0] - a[0], b[1] - a[1], b[2] - a[2]}};
    }
    md_util_min_image_vec3(row, m, cell);
    for (uint32_t j = 0; j < m; ++j) d[j] = sqrtf(row[j].x * row[j].x + row[j].y * row[j].y + row[j].z * row[j].z);
}

// The extreme pair of points P[n] and Q[m] under the minimum image, by the matrix's entries: the first pair in row
// order at the extreme of them wins ties, and the distance is its entry, the same bits distance() gives for that pair.
// Distinct leaves out j == skip and, given pp and qp, pairs of one particle. bi is -1 when no pair is left. row and d
// hold m each. One definition for a call and for a map's calls at once (vial_k_pair_batch), so both give the same
// bits; pair_extreme_fast finds the same pair faster.
static void pair_extreme(const float* P, uint32_t n, const float* Q, uint32_t m, const md_unitcell_t* cell, bool closest, bool distinct,
                         int64_t skip, const uint32_t* pp, const uint32_t* qp, vec3_t* row, float* d, float* best_out, int32_t* bi_out, int32_t* bj_out) {
    float best = 0.0f;
    int32_t bi = -1, bj = -1;
    for (uint32_t i = 0; i < n; ++i) {
        pair_row(P + 3 * i, Q, m, cell, row, d);
        for (uint32_t j = 0; j < m; ++j) {
            if (distinct && ((int64_t)j == skip || (pp && qp && pp[i] == qp[j]))) continue;
            if (bi < 0 || (closest ? d[j] < best : d[j] > best)) { best = d[j]; bi = (int32_t)i; bj = (int32_t)j; }
        }
    }
    *best_out = best;
    *bi_out = bi;
    *bj_out = bj;
}

// The largest magnitude of v[n], false if one is not finite. Eight at a time: positions of every particle of a
// trajectory pass through it each frame.
static bool finite_scale(const float* v, size_t n, float* scale) {
    const md_256 abs_mask = md_mm256_castsi256_ps(md_mm256_set1_epi32(0x7FFFFFFF));
    const md_256 fmax = md_mm256_set1_ps(FLT_MAX);
    md_256 s = md_mm256_set1_ps(*scale);
    int bad = 0;
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        const md_256 a = md_mm256_and_ps(md_mm256_loadu_ps(v + i), abs_mask);
        // A NaN or an infinity is not <= FLT_MAX
        bad |= md_mm256_movemask_ps(md_mm256_cmple_ps(a, fmax)) ^ 0xFF;
        s = md_mm256_max_ps(s, a);
    }
    float m = md_mm256_reduce_max_ps(s);
    for (; i < n; ++i) {
        const float a = fabsf(v[i]);
        bad |= !(a <= FLT_MAX);
        m = a > m ? a : m;
    }
    *scale = m;
    return !bad;
}

// pair_extreme without distinct, by mdlib's distance kernel, which skips what cannot be the extreme: the same pair
// and the same bits. The kernel squares and reduces in another order than the matrix's entries, so its distances may
// be a few ulps off them, within slack, a bound with ample room. The pair at the extreme of the entries is then within
// twice that of the kernel's extreme, so the kernel collects every pair within it, in the same pass, and the extreme
// of their entries is the matrix's. False where it does not apply: small inputs, which pair_extreme does as fast, a
// triclinic cell, whose minimum image the kernel takes otherwise than distance() does, positions that are not
// finite, and more near pairs than a tie of a few thousand.
static bool pair_extreme_fast(md_allocator_i* alloc, const float* P, uint32_t n, const float* Q, uint32_t m, const md_unitcell_t* cell,
                              bool closest, float* best_out, int32_t* bi_out, int32_t* bj_out) {
    if ((uint64_t)n * m < 64 || (md_unitcell_flags(cell) & MD_UNITCELL_TRICLINIC)) return false;
    float scale = 0.0f;
    if (!finite_scale(P, 3 * (size_t)n, &scale) || !finite_scale(Q, 3 * (size_t)m, &scale)) return false;
    float ext[3];
    md_unitcell_diag_extract_float(ext, cell);
    const float slack = (scale + ext[0] + ext[1] + ext[2]) * 1.0e-5f + 1.0e-4f;
    for (size_t cap = 64; cap <= 4096; cap *= 8) {
        int64_t* a = md_alloc(alloc, cap * sizeof(int64_t));
        int64_t* b = md_alloc(alloc, cap * sizeof(int64_t));
        const size_t cnt = closest ? md_util_min_distance_near(a, b, cap, 2.0f * slack, (const vec3_t*)P, n, (const vec3_t*)Q, m, cell)
                                   : md_util_max_distance_near(a, b, cap, 2.0f * slack, (const vec3_t*)P, n, (const vec3_t*)Q, m, cell);
        if (cnt > cap) continue;
        float best = 0.0f;
        int32_t bi = -1, bj = -1;
        for (size_t c = 0; c < cnt; ++c) {
            const int32_t i = (int32_t)a[c], j = (int32_t)b[c];
            const float e = pair_entry(P + 3 * (size_t)i, Q + 3 * (size_t)j, cell);
            // The first in row order among equal entries, as pair_extreme has it
            if (bi < 0 || (closest ? e < best : e > best) || (e == best && (i < bi || (i == bi && j < bj)))) { best = e; bi = i; bj = j; }
        }
        if (bi < 0) return false;   // cannot be: the kernel's own extreme is near itself
        *best_out = best;
        *bi_out = bi;
        *bj_out = bj;
        return true;
    }
    return false;
}

bool vial_k_pair_batch(vial_kctx_t* k, bool closest, const vial_sets_t* A, const vial_sets_t* B, const uint32_t* ia, const uint32_t* ib,
                       size_t count, vial_val_t* out) {
    const float* xyz = frame_xyz(k);
    if (!xyz) return false;
    const size_t na = k->in->frame->num_atoms;
    const md_unitcell_t cell = frame_cell(k);
    // Every group's positions once, packed, as vial_k_positions gives them
    const vial_sets_t* G[2] = {A, B};
    float* pos[2];
    uint32_t mx[2] = {0, 0};
    for (int g = 0; g < 2; ++g) {
        const vial_sets_t* s = G[g];
        const uint32_t total = s->offset[s->count];
        pos[g] = md_alloc(k->alloc, (total ? total : 1) * 3 * sizeof(float));
        for (uint32_t t = 0; t < total; ++t) {
            if (s->idx[t] >= na) return vial_rt_fail(k, "R0103", "particle %u is beyond the frame's %zu", s->idx[t], na);
            memcpy(pos[g] + 3 * (size_t)t, xyz + 3 * (size_t)s->idx[t], 3 * sizeof(float));
        }
        for (uint32_t c = 0; c < s->count; ++c) if (vial_sets_size(s, c) > mx[g]) mx[g] = vial_sets_size(s, c);
    }
    vec3_t* row = md_alloc(k->alloc, (mx[1] ? mx[1] : 1) * sizeof(vec3_t));
    float* d2 = md_alloc(k->alloc, (mx[1] ? mx[1] : 1) * sizeof(float));
    float* O = out->data;
    for (size_t c = 0; c < count; ++c) {
        if ((c & 255) == 0 && k->in->stop && md_atomic_load64(k->in->stop)) { k->in->interrupted = true; return false; }
        const uint32_t n = vial_sets_size(A, ia[c]), m = vial_sets_size(B, ib[c]);
        float best = 0.0f;
        int32_t bi = -1, bj = -1;
        if (n && m) {
            const float* P = pos[0] + 3 * (size_t)A->offset[ia[c]];
            const float* Q = pos[1] + 3 * (size_t)B->offset[ib[c]];
            if (k->in->plain || !pair_extreme_fast(k->alloc, P, n, Q, m, &cell, closest, &best, &bi, &bj)) {
                pair_extreme(P, n, Q, m, &cell, closest, false, -1, NULL, NULL, row, d2, &best, &bi, &bj);
            }
        }
        if (bi < 0) { O[c] = 0.0f; vial_val_set_missing(out, c, k->alloc); }
        else O[c] = best;
    }
    return true;
}

static bool k_pair_x(vial_kctx_t* k, const vial_val_t* pv, const vial_val_t* qv, bool closest, bool distinct, int64_t skip, vial_val_t* out) {
    vial_val_t p, q;
    if (!points(k, pv, &p) || !points(k, qv, &q)) return false;
    if (p.ndim != 2 || q.ndim != 2) return vial_rt_fail(k, "R0103", "closest_pair over arrays of groups is mapped by its def");
    const uint32_t n = p.dims[0], m = q.dims[0];
    const md_unitcell_t cell = frame_cell(k);
    vial_val_t* el = md_alloc(k->alloc, 3 * sizeof(vial_val_t));
    el[0] = vial_val_float(0.0f, k->alloc);
    el[1] = vial_val_int(0, k->alloc);
    el[2] = vial_val_int(0, k->alloc);
    if (!n || !m) {
        vial_val_set_missing(&el[0], 0, k->alloc);
        vial_val_set_missing(&el[1], 0, k->alloc);
        vial_val_set_missing(&el[2], 0, k->alloc);
    } else {
        const bool by_particle = distinct && p.origin && q.origin;
        vec3_t* row = md_alloc(k->alloc, m * sizeof(vec3_t));
        float* d2 = md_alloc(k->alloc, m * sizeof(float));
        float best = 0.0f;
        int32_t bi = -1, bj = -1;
        if (distinct || p.missing || q.missing || k->in->plain ||
            !pair_extreme_fast(k->alloc, (const float*)p.data, n, (const float*)q.data, m, &cell, closest, &best, &bi, &bj)) {
            pair_extreme((const float*)p.data, n, (const float*)q.data, m, &cell, closest, distinct, skip,
                by_particle ? vial_sets_items(p.origin, 0) : NULL, by_particle ? vial_sets_items(q.origin, 0) : NULL, row, d2, &best, &bi, &bj);
        }
        *(float*)el[0].data = best;
        *(int32_t*)el[1].data = bi;
        *(int32_t*)el[2].data = bj;
        if (bi < 0) for (int e = 0; e < 3; ++e) vial_val_set_missing(&el[e], 0, k->alloc);
    }
    *out = (vial_val_t){.dtype = VIAL_DT_TUPLE, .elems = el, .num_elems = 3};
    return true;
}

// ### PROPERTIES ###

static bool k_property(vial_kctx_t* k, const char* name, const vial_val_t* s, vial_val_t* out) {
    if (!s) return vial_rt_fail(k, "R0103", "a missing argument");
    const md_system_t* sys = k->in->msys;
    if (!sys) return vial_rt_fail(k, "R0103", "%s needs a system", name);
    if (s->dtype != VIAL_DT_SET) {
        if (s->origin) return k_property(k, name, &(vial_val_t){.dtype = VIAL_DT_SET, .sets = s->origin}, out);
        return vial_rt_fail(k, "R0103", "%s of what is not a set", name);
    }
    if (s->ndim) {
        // An array of sets of one particle: one value per set
        const size_t cnt = vial_val_count(s);
        uint32_t* one = md_alloc(k->alloc, (cnt ? cnt : 1) * sizeof(uint32_t));
        for (size_t i = 0; i < cnt; ++i) {
            if (vial_sets_size(s->sets, (uint32_t)i) != 1) return vial_rt_fail(k, "R0103", "%s of sets of several particles", name);
            one[i] = vial_sets_items(s->sets, (uint32_t)i)[0];
        }
        vial_val_t flat;
        if (!k_property(k, name, &(vial_val_t){.dtype = VIAL_DT_SET, .sets = vial_sets_one(one, (uint32_t)cnt, k->alloc)}, &flat)) return false;
        flat.ndim = s->ndim;
        memcpy(flat.dims, s->dims, sizeof(s->dims));
        *out = flat;
        return true;
    }
    const uint32_t n = vial_sets_size(s->sets, 0);
    const uint32_t* idx = vial_sets_items(s->sets, 0);
    const bool str = !strcmp(name, "name") || !strcmp(name, "element");
    vial_val_t o; if (!vial_k_make(k, str ? VIAL_DT_STR : VIAL_DT_FLOAT, 1, &n, &o)) return false;
    const md_attribute_t* charges = !strcmp(name, "charge") ? md_attributes_find(&sys->attributes, STR_LIT("atom/charge")) : NULL;
    const float* cq = charges ? md_attribute_view(charges, MD_ATTRIBUTE_TYPE_F32, 1, 1) : NULL;
    const int which = !strcmp(name, "mass") ? 0 : !strcmp(name, "charge") ? 1 : !strcmp(name, "name") ? 2 : !strcmp(name, "element") ? 3 : -1;
    if (which < 0) return vial_rt_fail(k, "R0103", "%s is not evaluated yet", name);
    const float* ms = which == 0 ? k_masses(k) : NULL;
    for (uint32_t i = 0; i < n; ++i) {
        const size_t a = idx[i];
        switch (which) {
        case 0:  ((float*)o.data)[i] = k_mass(ms, sys, a); break;
        case 1:  ((float*)o.data)[i] = cq ? cq[a] : (float)md_atom_formal_charge(&sys->atom, a); break;
        case 2:  ((str_t*)o.data)[i] = md_atom_name(&sys->atom, a); break;
        default: ((str_t*)o.data)[i] = md_util_element_symbol(md_atom_atomic_number(&sys->atom, a)); break;
        }
    }
    o.origin = s->sets;
    *out = o;
    return true;
}

// ### REDUCTIONS ###

// Which dimensions a reduction takes: axis= as an int, 'frame', a tuple of them, or none for every dimension but the
// frame. On a series (VIAL_VAL_FRAMES) the frame is dimension 0 and the value's own dimensions follow it.
static bool axes_of(vial_kctx_t* k, const vial_val_t* axis, const vial_val_t* a, bool* red) {
    const uint8_t ndim = a->ndim;
    const bool frames = a->flags & VIAL_VAL_FRAMES;
    const uint8_t first = frames ? 1 : 0;
    for (uint8_t d = 0; d < ndim; ++d) red[d] = (axis == NULL || axis->dtype == VIAL_DT_NONE) && d >= first;
    if (!axis || axis->dtype == VIAL_DT_NONE) return true;
    const vial_val_t* items = axis->dtype == VIAL_DT_TUPLE ? axis->elems : axis;
    const uint32_t n = axis->dtype == VIAL_DT_TUPLE ? axis->num_elems : 1;
    for (uint32_t i = 0; i < n; ++i) {
        if (items[i].dtype == VIAL_DT_STR) {
            if (!frames) return vial_rt_fail(k, "R0103", "axis='frame' reduces a series of frames");
            red[0] = true;
            continue;
        }
        if (items[i].dtype != VIAL_DT_INT) return vial_rt_fail(k, "R0103", "axis= takes ints or 'frame'");
        int64_t ax = *(const int32_t*)items[i].data;
        const int64_t own = ndim - first;
        if (ax < 0) ax += own;
        if (ax < 0 || ax >= own) return vial_rt_fail(k, "R0103", "axis %d of %lld dimensions", *(const int32_t*)items[i].data, (long long)own);
        red[ax + first] = true;
    }
    return true;
}

typedef enum { R_SUM, R_MEAN, R_MIN, R_MAX, R_STD, R_VAR, R_ANY, R_ALL, R_ARGMIN, R_ARGMAX } red_t;

// One reduction over a strided run of n values. Sums go to 8 lanes in order, then a pairwise tree, in 64 bits: the
// same order on every machine, whatever its vector width.
typedef struct racc_t {
    double lane[8], lane_sin[8], lane_cos[8], lane_sq[8];
    size_t count;
    double best;
    int64_t best_at;
    bool any_nan;
} racc_t;

static double tree8(const double* l) { return ((l[0] + l[1]) + (l[2] + l[3])) + ((l[4] + l[5]) + (l[6] + l[7])); }

static bool k_reduce(vial_kctx_t* k, red_t red, const vial_val_t* a, const vial_val_t* axis, vial_val_t* out) {
    if (!a) return vial_rt_fail(k, "R0103", "a missing argument");
    const vial_type_t* rt = k->type;
    if (a->dtype == VIAL_DT_SET) {
        // any and all of sets: their union and intersection
        if (red != R_ANY && red != R_ALL) return vial_rt_fail(k, "R0103", "a reduction of sets other than any and all");
        const uint32_t n = a->sets->count;
        vial_sets_t* acc = NULL;
        for (uint32_t i = 0; i < n; ++i) {
            if (!acc) { acc = vial_sets_pick(a->sets, &i, 1, k->alloc); continue; }
            acc = red == R_ANY ? vial_sets_or(acc, 0, a->sets, i, k->alloc) : vial_sets_and(acc, 0, a->sets, i, k->alloc);
        }
        if (!acc) acc = vial_sets_one(NULL, 0, k->alloc);
        acc->single = true;
        *out = vial_val_sets(acc);
        return true;
    }
    if (!is_num(a->dtype)) return vial_rt_fail(k, "R0103", "a reduction of what is not a number");
    bool rd[VIAL_MAX_DIMS];
    if (!axes_of(k, axis, a, rd)) return false;
    // Output dimensions: those not reduced; the reduced ones are walked for each output element
    uint32_t odims[VIAL_MAX_DIMS], rdims[VIAL_MAX_DIMS];
    size_t ostr[VIAL_MAX_DIMS], rstr[VIAL_MAX_DIMS], st[VIAL_MAX_DIMS];
    uint8_t ond = 0, rnd = 0;
    {
        size_t s = 1;
        for (int d = a->ndim - 1; d >= 0; --d) { st[d] = s; s *= a->dims[d]; }
    }
    for (uint8_t d = 0; d < a->ndim; ++d) {
        if (rd[d]) { rdims[rnd] = a->dims[d]; rstr[rnd++] = st[d]; }
        else { odims[ond] = a->dims[d]; ostr[ond++] = st[d]; }
    }
    size_t on = 1, rn = 1;
    for (uint8_t d = 0; d < ond; ++d) on *= odims[d];
    for (uint8_t d = 0; d < rnd; ++d) rn *= rdims[d];
    const bool periodic = (red == R_MEAN || red == R_STD || red == R_VAR) && k->op && k->op->num_args &&
                          (k->in->ir->values[k->in->ir->args[k->op->arg_beg]].type.flags & VIAL_TYPE_PERIODIC);
    const vial_type_t* at = k->op && k->op->num_args ? &k->in->ir->values[k->in->ir->args[k->op->arg_beg]].type : NULL;
    const double pw = periodic ? (double)(at->period_end - at->period_beg) : 0.0;
    const uint8_t odt = rt->dtype;
    vial_val_t o; if (!vial_k_make(k, odt, ond, odims, &o)) return false;
    for (size_t oi = 0; oi < on; ++oi) {
        size_t base = 0, rem = oi;
        for (int d = ond - 1; d >= 0; --d) { base += (rem % odims[d]) * ostr[d]; rem /= odims[d]; }
        racc_t acc;
        memset(&acc, 0, sizeof(acc));
        acc.best_at = -1;
        bool bool_acc = red == R_ALL;
        size_t j = 0;
        for (size_t ri = 0; ri < rn; ++ri) {
            size_t off = base, r2 = ri;
            for (int d = rnd - 1; d >= 0; --d) { off += (r2 % rdims[d]) * rstr[d]; r2 /= rdims[d]; }
            if (vial_val_is_missing(a, off)) continue;     // reductions skip missing entries
            const double v = num_at(a, off);
            const int lane = (int)(j & 7);
            j += 1;
            acc.count += 1;
            switch (red) {
            case R_SUM: case R_MEAN: case R_STD: case R_VAR:
                if (periodic) {
                    const double rad = (v - at->period_beg) / pw * 2.0 * 3.141592653589793;
                    acc.lane_sin[lane] += sin(rad);
                    acc.lane_cos[lane] += cos(rad);
                } else {
                    acc.lane[lane] += v;
                    acc.lane_sq[lane] += v * v;
                }
                break;
            case R_MIN: case R_ARGMIN:
                if (v != v) acc.any_nan = true;
                if (acc.best_at < 0 || v < acc.best) { acc.best = v; acc.best_at = (int64_t)ri; }
                break;
            case R_MAX: case R_ARGMAX:
                if (v != v) acc.any_nan = true;
                if (acc.best_at < 0 || v > acc.best) { acc.best = v; acc.best_at = (int64_t)ri; }
                break;
            case R_ANY: bool_acc = bool_acc || v != 0.0; break;
            case R_ALL: bool_acc = bool_acc && v != 0.0; break;
            }
        }
        double r = 0.0;
        bool miss = false;
        switch (red) {
        case R_SUM: r = tree8(acc.lane); break;
        case R_MEAN: case R_STD: case R_VAR:
            if (!acc.count) { miss = true; break; }
            if (periodic) {
                const double sx = tree8(acc.lane_sin), cx = tree8(acc.lane_cos);
                const double R = sqrt(sx * sx + cx * cx) / (double)acc.count;
                if (R <= 0.0) { miss = true; break; }
                if (red == R_MEAN) {
                    double m = atan2(sx, cx) / (2.0 * 3.141592653589793) * pw + at->period_beg;
                    while (m <= at->period_beg) m += pw;
                    while (m > at->period_end) m -= pw;
                    r = m;
                } else {
                    const double sd = sqrt(-2.0 * log(R > 1.0 ? 1.0 : R)) * pw / (2.0 * 3.141592653589793);
                    r = red == R_STD ? sd : sd * sd;
                }
            } else {
                const double mean = tree8(acc.lane) / (double)acc.count;
                if (red == R_MEAN) r = mean;
                else {
                    double var = tree8(acc.lane_sq) / (double)acc.count - mean * mean;
                    if (var < 0.0) var = 0.0;
                    r = red == R_VAR ? var : sqrt(var);
                }
            }
            break;
        case R_MIN: case R_MAX:
            if (!acc.count) miss = true;
            r = acc.any_nan ? NAN : acc.best;
            break;
        case R_ARGMIN: case R_ARGMAX:
            if (!acc.count) miss = true;
            r = (double)acc.best_at;
            break;
        case R_ANY: case R_ALL: r = bool_acc; break;
        }
        switch (odt) {
        case VIAL_DT_BOOL:    ((uint8_t*)o.data)[oi] = r != 0.0; break;
        case VIAL_DT_INT:
            if (r > INT32_MAX || r < INT32_MIN) miss = true;     // R0102
            else ((int32_t*)o.data)[oi] = (int32_t)r;
            break;
        case VIAL_DT_FLOAT:   ((float*)o.data)[oi] = (float)r; break;
        case VIAL_DT_FLOAT64: ((double*)o.data)[oi] = r; break;
        default: return vial_rt_fail(k, "R0103", "a reduction to what is not a number");
        }
        if (miss) vial_val_set_missing(&o, oi, k->alloc);
    }
    *out = o;
    return true;
}

// ### REDUCTIONS OVER EVERY PAIR ###
// VIAL_IR_PAIRS: a line of the matrix at a time, a row or a column, computed as displacement and norm compute an entry
// (b - a, mdlib's minimum image over the line, the squares summed in order and sqrtf) and reduced by k_reduce itself,
// so every result has the matrix's bits. Over every pair, each row is reduced and then the rows' results, which min,
// max, a count and any allow exactly. Memory is a line, not n * m.

// A count or any of d < r over the pairs the spatial search finds within a hair past r, each tested with the matrix's
// own entry: the same pairs pass as in the matrix, without computing the others.
typedef struct pairs_search_t {
    const float*        px;
    const float*        py;
    const md_unitcell_t* cell;
    uint32_t            test;
    double              r;
    bool                rows;
    int32_t*            counts;     // per row or column
} pairs_search_t;

static void pairs_search_cb(const uint32_t* i_idx, const uint32_t* j_idx, const float* d2, size_t n, void* user) {
    (void)d2;
    pairs_search_t* s = user;
    for (size_t q = 0; q < n; ++q) {
        const uint32_t i = i_idx[q], j = j_idx[q];
        const float* A = s->px + 3 * (size_t)i;
        const float* B = s->py + 3 * (size_t)j;
        vec3_t v = {{B[0] - A[0], B[1] - A[1], B[2] - A[2]}};
        md_util_min_image_vec3(&v, 1, s->cell);
        float e = 0.0f;
        e += v.x * v.x;
        e += v.y * v.y;
        e += v.z * v.z;
        const float d = sqrtf(e);
        const bool pass = s->test == VIAL_CMP_LT ? (double)d < s->r : (double)d <= s->r;
        if (pass) s->counts[s->rows ? i : j] += 1;
    }
}

// Whether the search reports every pair within reach once, as one image: no triclinic cell, and reach short of half
// of every orthorhombic period
static bool search_ok(const md_unitcell_t* cell, double reach) {
    if (cell->flags & MD_UNITCELL_TRICLINIC) return false;
    if (!(cell->flags & MD_UNITCELL_PBC_ALL)) return true;
    float ext[3];
    md_unitcell_diag_extract_float(ext, cell);
    const uint32_t pbc[3] = {MD_UNITCELL_PBC_X, MD_UNITCELL_PBC_Y, MD_UNITCELL_PBC_Z};
    for (int c = 0; c < 3; ++c) if ((cell->flags & pbc[c]) && !(reach < 0.5 * ext[c])) return false;
    return true;
}

static bool all_finite(const float* v, size_t n) {
    for (size_t i = 0; i < n; ++i) if (!isfinite(v[i])) return false;
    return true;
}

// The nearest of each line among the pairs the search finds, each with the matrix's own entry; ties to the first
typedef struct pairs_nearest_t {
    const float*        px;
    const float*        py;
    const md_unitcell_t* cell;
    bool                rows;
    float*              best;
    int32_t*            at;
} pairs_nearest_t;

static void pairs_nearest_cb(const uint32_t* i_idx, const uint32_t* j_idx, const float* d2, size_t n, void* user) {
    (void)d2;
    pairs_nearest_t* s = user;
    for (size_t q = 0; q < n; ++q) {
        const uint32_t i = i_idx[q], j = j_idx[q];
        const float* A = s->px + 3 * (size_t)i;
        const float* B = s->py + 3 * (size_t)j;
        vec3_t v = {{B[0] - A[0], B[1] - A[1], B[2] - A[2]}};
        md_util_min_image_vec3(&v, 1, s->cell);
        float e = 0.0f;
        e += v.x * v.x;
        e += v.y * v.y;
        e += v.z * v.z;
        const float d = sqrtf(e);
        const uint32_t L = s->rows ? i : j, other = s->rows ? j : i;
        if (s->at[L] < 0 || d < s->best[L] || (d == s->best[L] && (int32_t)other < s->at[L])) { s->best[L] = d; s->at[L] = (int32_t)other; }
    }
}

// within(): every particle against groups of eight points of `of`, rather than a spatial search over every particle,
// up to this many tests of eight particles against a group, about a millisecond
#define VIAL_WITHIN_BRUTE 262144

// How far the search for nearest pairs looks: a line whose nearest is past it is computed whole
#define VIAL_PAIRS_NEAR 8.0

bool vial_k_pairs(vial_kctx_t* k, uint32_t red, uint32_t axis, uint32_t test, const vial_val_t* x, const vial_val_t* y, const vial_val_t* r,
                  vial_val_t* out) {
    if (!x || !y || (test && !r)) return vial_rt_fail(k, "R0103", "a missing argument");
    vial_val_t px, py;
    if (!point(k, x, &px) || !point(k, y, &py)) return false;
    const size_t n = vial_val_count(&px) / 3, m = vial_val_count(&py) / 3;
    const red_t rd = red == VIAL_PAIRS_MIN ? R_MIN : red == VIAL_PAIRS_MAX ? R_MAX : red == VIAL_PAIRS_ARGMIN ? R_ARGMIN :
                     red == VIAL_PAIRS_ARGMAX ? R_ARGMAX : red == VIAL_PAIRS_COUNT ? R_SUM : R_ANY;
    const uint8_t odt = k->type->dtype;
    const md_unitcell_t cell = frame_cell(k);
    const double rv = test ? num_at(r, 0) : 0.0;
    const bool r_missing = test && vial_val_is_missing(r, 0);
    // The lines: rows of m entries, or columns of n
    const bool rows = axis != 0;

    // A count or any within a cutoff: the spatial search, where every pair it reports is one image, the cutoff short
    // of half an orthorhombic period, and no position is missing
    if (test && (red == VIAL_PAIRS_COUNT || red == VIAL_PAIRS_ANY) && !px.missing && !py.missing && !r_missing && rv >= 0.0 && rv < 1e6 && n && m) {
        const double reach = nextafter(rv, 1e30) * 1.00001 + 1e-4;
        if (search_ok(&cell, reach)) {
            int32_t* counts = md_alloc(k->alloc, (rows ? n : m) * sizeof(int32_t));
            memset(counts, 0, (rows ? n : m) * sizeof(int32_t));
            pairs_search_t ps = {.px = px.data, .py = py.data, .cell = &cell, .test = test, .r = rv, .rows = rows, .counts = counts};
            md_coord_stream_t ys = md_coord_stream_from_aos(py.data, 3 * sizeof(float), NULL, m);
            md_coord_stream_t xs = md_coord_stream_from_aos(px.data, 3 * sizeof(float), NULL, n);
            md_spatial_acc_t acc = {.alloc = k->alloc};
            k_spatial_init(k, &acc, &ys, reach, &cell);
            md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &xs, reach, pairs_search_cb, &ps, 0);
            md_spatial_acc_free(&acc);
            const size_t nl = rows ? n : m;
            if (axis == VIAL_PAIRS_ALL) {
                int64_t total = 0;
                for (size_t L = 0; L < nl; ++L) total += counts[L];
                vial_val_t o;
                if (!vial_k_make(k, odt, 0, NULL, &o)) return false;
                if (red == VIAL_PAIRS_ANY) *(uint8_t*)o.data = total > 0;
                else if (total > INT32_MAX) vial_val_set_missing(&o, 0, k->alloc);     // R0102, as the sum gives
                else *(int32_t*)o.data = (int32_t)total;
                *out = o;
                return true;
            }
            vial_val_t o;
            const uint32_t od[1] = {(uint32_t)nl};
            if (!vial_k_make(k, odt, 1, od, &o)) return false;
            for (size_t L = 0; L < nl; ++L) {
                if (red == VIAL_PAIRS_ANY) ((uint8_t*)o.data)[L] = counts[L] > 0;
                else ((int32_t*)o.data)[L] = counts[L];
            }
            *out = o;
            return true;
        }
    }
    // The minimum or maximum over every pair: the closest or farthest pair's entry, by the kernels
    if (!test && axis == VIAL_PAIRS_ALL && (red == VIAL_PAIRS_MIN || red == VIAL_PAIRS_MAX) && odt == VIAL_DT_FLOAT && !px.missing && !py.missing && n && m && !k->in->plain) {
        float best;
        int32_t bi, bj;
        if (pair_extreme_fast(k->alloc, px.data, (uint32_t)n, py.data, (uint32_t)m, &cell, red == VIAL_PAIRS_MIN, &best, &bi, &bj)) {
            vial_val_t o;
            if (!vial_k_make(k, odt, 0, NULL, &o)) return false;
            *(float*)o.data = best;
            *out = o;
            return true;
        }
    }
    const size_t lines = rows ? n : m, len = rows ? m : n;
    vec3_t* buf = md_alloc(k->alloc, (len ? len : 1) * sizeof(vec3_t));
    vial_val_t line = {.dtype = test ? VIAL_DT_BOOL : VIAL_DT_FLOAT, .ndim = 1, .dims = {(uint32_t)len}};
    line.data = md_alloc(k->alloc, (len ? len : 1) * sizeof(float));
    const size_t mask_bytes = ((len ? len : 1) + 63) / 64 * sizeof(uint64_t);
    uint64_t* mask = md_alloc(k->alloc, mask_bytes);
    // One line's reduction to a scalar of the result's dtype, or of the partial's over every pair
    const vial_type_t lt = vial_type_scalar((vial_dtype_t)odt);
    vial_kctx_t lk = *k;
    lk.type = &lt;
    vial_val_t o;
    const uint32_t odims[1] = {(uint32_t)lines};
    if (!vial_k_make(k, odt, 1, odims, &o)) return false;
    const size_t es = vial_dtype_size(odt);

    // The nearest of each line: every pair within VIAL_PAIRS_NEAR is found by the search, so a line whose nearest found
    // is that close has it exactly, ties to the first as argmin has them; the other lines are computed whole below
    uint8_t* resolved = NULL;
    const double near_reach = VIAL_PAIRS_NEAR * 1.00001 + 1e-4;
    if (!test && (red == VIAL_PAIRS_MIN || red == VIAL_PAIRS_ARGMIN) && !px.missing && !py.missing && n && m && search_ok(&cell, near_reach) &&
        all_finite(px.data, 3 * n) && all_finite(py.data, 3 * m)) {
        float* best = md_alloc(k->alloc, lines * sizeof(float));
        int32_t* at = md_alloc(k->alloc, lines * sizeof(int32_t));
        for (size_t L = 0; L < lines; ++L) { best[L] = 0.0f; at[L] = -1; }
        pairs_nearest_t ns = {.px = px.data, .py = py.data, .cell = &cell, .rows = rows, .best = best, .at = at};
        md_coord_stream_t ys = md_coord_stream_from_aos(py.data, 3 * sizeof(float), NULL, m);
        md_coord_stream_t xs = md_coord_stream_from_aos(px.data, 3 * sizeof(float), NULL, n);
        md_spatial_acc_t acc = {.alloc = k->alloc};
        k_spatial_init(k, &acc, &ys, near_reach, &cell);
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &xs, near_reach, pairs_nearest_cb, &ns, 0);
        md_spatial_acc_free(&acc);
        resolved = md_alloc(k->alloc, lines);
        for (size_t L = 0; L < lines; ++L) {
            resolved[L] = at[L] >= 0 && (double)best[L] <= VIAL_PAIRS_NEAR;
            if (!resolved[L]) continue;
            if (red == VIAL_PAIRS_MIN) ((float*)o.data)[L] = best[L];
            else ((int32_t*)o.data)[L] = at[L];
        }
    }
    for (size_t L = 0; L < lines; ++L) {
        if ((L & 63) == 0 && k->in->stop && md_atomic_load64(k->in->stop)) { k->in->interrupted = true; return false; }
        if (resolved && resolved[L]) continue;
        bool any_missing = r_missing;
        for (size_t e = 0; e < len; ++e) {
            const size_t i = rows ? L : e, j = rows ? e : L;
            const float* A = (const float*)px.data + 3 * i;
            const float* B = (const float*)py.data + 3 * j;
            buf[e] = (vec3_t){{B[0] - A[0], B[1] - A[1], B[2] - A[2]}};
        }
        md_util_min_image_vec3(buf, len, &cell);
        memset(mask, 0, mask_bytes);
        for (size_t e = 0; e < len; ++e) {
            const size_t i = rows ? L : e, j = rows ? e : L;
            float s = 0.0f;
            s += buf[e].x * buf[e].x;
            s += buf[e].y * buf[e].y;
            s += buf[e].z * buf[e].z;
            const float d = sqrtf(s);
            bool miss = r_missing;
            for (int c = 0; c < 3 && !miss; ++c) miss = vial_val_is_missing(&px, 3 * i + c) || vial_val_is_missing(&py, 3 * j + c);
            if (miss) { mask[e >> 6] |= 1ull << (e & 63); any_missing = true; }
            if (test) ((uint8_t*)line.data)[e] = test == VIAL_CMP_LT ? (double)d < rv : (double)d <= rv;
            else ((float*)line.data)[e] = d;
        }
        line.missing = any_missing ? mask : NULL;
        vial_val_t res;
        if (!k_reduce(&lk, rd, &line, NULL, &res)) return false;
        memcpy((uint8_t*)o.data + L * es, res.data, es);
        if (vial_val_is_missing(&res, 0)) vial_val_set_missing(&o, L, k->alloc);
    }
    if (axis != VIAL_PAIRS_ALL) { *out = o; return true; }
    // Over every pair: the rows' results reduced, as the matrix's entries are; a count sums the rows' counts
    return k_reduce(k, rd, &o, NULL, out);
}

// ### ELEMENTWISE MATH ###

static bool k_math1(vial_kctx_t* k, const char* name, const vial_val_t* a, vial_val_t* out) {
    if (!a) return vial_rt_fail(k, "R0103", "a missing argument");
    const uint8_t odt = k->type->dtype;
    const vial_type_t* at = k->op && k->op->num_args ? &k->in->ir->values[k->in->ir->args[k->op->arg_beg]].type : NULL;
    // An angle in degrees is turned to radians for sin, cos and tan
    double scale = 1.0;
    if (at && at->unit_kind == VIAL_UNIT_KNOWN && (!strcmp(name, "sin") || !strcmp(name, "cos") || !strcmp(name, "tan"))) {
        double f;
        if (md_unit_conversion_factor(&f, at->unit, md_unit_radian())) scale = f;
    }
    const size_t n = vial_val_count(a);
    if (odt == VIAL_DT_INT) {
        const vial_val_t x = vial_val_as(a, VIAL_DT_INT, k->alloc);
        vial_val_t o; if (!vial_k_make(k, VIAL_DT_INT, x.ndim, x.dims, &o)) return false;
        const bool is_abs = !strcmp(name, "abs");
        for (size_t i = 0; i < n; ++i) {
            const int32_t v = ((const int32_t*)x.data)[i];
            if (is_abs) { if (v == INT32_MIN) vial_val_set_missing(&o, i, k->alloc); else ((int32_t*)o.data)[i] = v < 0 ? -v : v; }
            else ((int32_t*)o.data)[i] = v;    // floor, ceil and round of ints
        }
        o.missing = o.missing ? o.missing : x.missing;
        *out = o;
        return true;
    }
    // In double, rounded once to the result's dtype: a float's sin is the float nearest the sine of its value, not
    // whatever sinf makes of it, and an angle in degrees is turned to radians before any rounding
    enum { M_ABS, M_FLOOR, M_CEIL, M_ROUND, M_SQRT, M_EXP, M_LOG, M_LOG10, M_SIN, M_COS, M_TAN, M_ASIN, M_ACOS, M_ATAN, M_COUNT };
    static const char* names[M_COUNT] = {"abs", "floor", "ceil", "round", "sqrt", "exp", "log", "log10", "sin", "cos", "tan", "asin", "acos", "atan"};
    int fn = -1;
    for (int f = 0; f < M_COUNT; ++f) if (!strcmp(name, names[f])) fn = f;
    if (fn < 0) return vial_rt_fail(k, "R0103", "%s is not evaluated yet", name);
    const vial_val_t x = vial_val_as(a, odt, k->alloc);
    vial_val_t o; if (!vial_k_make(k, odt, x.ndim, x.dims, &o)) return false;
    for (size_t i = 0; i < n; ++i) {
        const double v = num_at(&x, i) * scale;
        double r;
        switch (fn) {
        case M_ABS:   r = fabs(v); break;
        case M_FLOOR: r = floor(v); break;
        case M_CEIL:  r = ceil(v); break;
        case M_ROUND: r = nearbyint(v); break;
        case M_SQRT:  r = sqrt(v); break;
        case M_EXP:   r = exp(v); break;
        case M_LOG:   r = log(v); break;
        case M_LOG10: r = log10(v); break;
        case M_SIN:   r = sin(v); break;
        case M_COS:   r = cos(v); break;
        case M_TAN:   r = tan(v); break;
        case M_ASIN:  r = asin(v) * DEG_PER_RAD; break;
        case M_ACOS:  r = acos(v) * DEG_PER_RAD; break;
        default:      r = atan(v) * DEG_PER_RAD; break;
        }
        if (odt == VIAL_DT_FLOAT) ((float*)o.data)[i] = (float)r;
        else ((double*)o.data)[i] = r;
    }
    o.missing = x.missing;
    o.origin = a->origin;
    *out = o;
    return true;
}

// Elementwise over broadcast arguments, in the result's dtype: where, clip, atan2, min and max of two
static bool k_elementwise(vial_kctx_t* k, const char* name, const vial_val_t* const* xs, uint32_t nx, vial_val_t* out) {
    const uint8_t odt = k->type->dtype;
    uint8_t nd = 0;
    for (uint32_t i = 0; i < nx; ++i) if (xs[i]->ndim > nd) nd = xs[i]->ndim;
    uint32_t dims[VIAL_MAX_DIMS];
    for (uint8_t d = 0; d < nd; ++d) dims[d] = 1;
    for (uint32_t i = 0; i < nx; ++i) {
        for (uint8_t d = 0; d < xs[i]->ndim; ++d) {
            uint32_t* o = &dims[nd - xs[i]->ndim + d];
            if (*o == 1) *o = xs[i]->dims[d];
            else if (xs[i]->dims[d] != 1 && xs[i]->dims[d] != *o) return vial_rt_fail(k, "R0103", "extents %u and %u do not broadcast", *o, xs[i]->dims[d]);
        }
    }
    size_t n = 1;
    for (uint8_t d = 0; d < nd; ++d) n *= dims[d];
    vial_val_t o; if (!vial_k_make(k, odt, nd, dims, &o)) return false;
    for (size_t i = 0; i < n; ++i) {
        double v[4] = {0};
        bool miss = false;
        for (uint32_t q = 0; q < nx; ++q) {
            size_t rem = i, off = 0, s = 1;
            for (int d = nd - 1; d >= 0; --d) {
                const uint32_t x = (uint32_t)(rem % dims[d]);
                rem /= dims[d];
                const int dd = d - (nd - xs[q]->ndim);
                if (dd >= 0) { if (xs[q]->dims[dd] != 1) off += x * s; s *= xs[q]->dims[dd]; }
            }
            v[q] = num_at(xs[q], off);
            miss = miss || vial_val_is_missing(xs[q], off);
        }
        double r;
        if (!strcmp(name, "where")) {
            // Only the chosen side's missing counts
            r = v[0] != 0.0 ? v[1] : v[2];
        } else if (!strcmp(name, "clip")) {
            r = v[0] < v[1] ? v[1] : v[0] > v[2] ? v[2] : v[0];
        } else if (!strcmp(name, "atan2")) {
            r = odt == VIAL_DT_FLOAT ? (double)(atan2f((float)v[0], (float)v[1]) * (float)DEG_PER_RAD) : atan2(v[0], v[1]) * DEG_PER_RAD;
        } else if (!strcmp(name, "min")) {
            r = v[0] != v[0] || v[1] != v[1] ? NAN : v[0] < v[1] ? v[0] : v[1];
        } else if (!strcmp(name, "max")) {
            r = v[0] != v[0] || v[1] != v[1] ? NAN : v[0] > v[1] ? v[0] : v[1];
        } else {
            return vial_rt_fail(k, "R0103", "%s is not evaluated yet", name);
        }
        switch (odt) {
        case VIAL_DT_BOOL:    ((uint8_t*)o.data)[i] = r != 0.0; break;
        case VIAL_DT_INT:     ((int32_t*)o.data)[i] = (int32_t)r; break;
        case VIAL_DT_FLOAT:   ((float*)o.data)[i] = (float)r; break;
        case VIAL_DT_FLOAT64: ((double*)o.data)[i] = r; break;
        }
        if (miss) vial_val_set_missing(&o, i, k->alloc);
    }
    *out = o;
    return true;
}

// norm and dot over the last dimension, in float with a fixed order
static bool k_norm(vial_kctx_t* k, const vial_val_t* a, bool normalize, vial_val_t* out) {
    if (!a) return vial_rt_fail(k, "R0103", "a missing argument");
    const vial_val_t x = vial_val_as(a, VIAL_DT_FLOAT, k->alloc);
    if (!x.ndim) return vial_rt_fail(k, "R0103", "norm of a scalar");
    const uint32_t m = x.dims[x.ndim - 1];
    const size_t g = m ? vial_val_count(&x) / m : 0;
    vial_val_t o;
    if (!vial_k_make(k, VIAL_DT_FLOAT, normalize ? x.ndim : x.ndim - 1, x.dims, &o)) return false;
    for (size_t i = 0; i < g; ++i) {
        const float* v = (const float*)x.data + i * m;
        float s = 0.0f;
        bool miss = false;
        for (uint32_t c = 0; c < m; ++c) { s += v[c] * v[c]; miss = miss || vial_val_is_missing(&x, i * m + c); }
        const float l = sqrtf(s);
        if (normalize) {
            for (uint32_t c = 0; c < m; ++c) {
                ((float*)o.data)[i * m + c] = v[c] / l;
                if (miss) vial_val_set_missing(&o, i * m + c, k->alloc);
            }
        } else {
            ((float*)o.data)[i] = l;
            if (miss) vial_val_set_missing(&o, i, k->alloc);
        }
    }
    if (!normalize) o.origin = a->origin;
    *out = o;
    return true;
}

static bool k_dot_cross(vial_kctx_t* k, const vial_val_t* a, const vial_val_t* b, bool cross, vial_val_t* out) {
    if (!a || !b) return vial_rt_fail(k, "R0103", "a missing argument");
    const vial_val_t x = vial_val_as(a, VIAL_DT_FLOAT, k->alloc), y = vial_val_as(b, VIAL_DT_FLOAT, k->alloc);
    if (!x.ndim || !y.ndim || x.dims[x.ndim - 1] != y.dims[y.ndim - 1]) return vial_rt_fail(k, "R0103", "last dimensions differ");
    const uint32_t m = x.dims[x.ndim - 1];
    if (cross && m != 3) return vial_rt_fail(k, "R0103", "cross takes a last dimension of 3");
    // Rows broadcast
    vial_val_t lx = x, ly = y;
    lx.ndim -= 1; ly.ndim -= 1;
    uint8_t nd = lx.ndim > ly.ndim ? lx.ndim : ly.ndim;
    uint32_t dims[VIAL_MAX_DIMS + 1];
    for (uint8_t d = 0; d < nd; ++d) dims[d] = 1;
    for (int w = 0; w < 2; ++w) {
        const vial_val_t* v = w ? &ly : &lx;
        for (uint8_t d = 0; d < v->ndim; ++d) {
            uint32_t* o = &dims[nd - v->ndim + d];
            if (*o == 1) *o = v->dims[d];
            else if (v->dims[d] != 1 && v->dims[d] != *o) return vial_rt_fail(k, "R0103", "rows do not broadcast");
        }
    }
    size_t n = 1;
    for (uint8_t d = 0; d < nd; ++d) n *= dims[d];
    if (cross) dims[nd] = 3;
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, cross ? nd + 1 : nd, dims, &o)) return false;
    for (size_t i = 0; i < n; ++i) {
        size_t rem = i, ia = 0, ib = 0, sa = 1, sb = 1;
        for (int d = nd - 1; d >= 0; --d) {
            const uint32_t q = (uint32_t)(rem % dims[d]);
            rem /= dims[d];
            const int da = d - (nd - lx.ndim), db = d - (nd - ly.ndim);
            if (da >= 0) { if (lx.dims[da] != 1) ia += q * sa; sa *= lx.dims[da]; }
            if (db >= 0) { if (ly.dims[db] != 1) ib += q * sb; sb *= ly.dims[db]; }
        }
        const float* u = (const float*)x.data + ia * m;
        const float* v = (const float*)y.data + ib * m;
        if (cross) {
            float* r = (float*)o.data + 3 * i;
            r[0] = u[1] * v[2] - u[2] * v[1];
            r[1] = u[2] * v[0] - u[0] * v[2];
            r[2] = u[0] * v[1] - u[1] * v[0];
        } else {
            float s = 0.0f;
            for (uint32_t c = 0; c < m; ++c) s += u[c] * v[c];
            ((float*)o.data)[i] = s;
        }
    }
    *out = o;
    return true;
}

// ### ARRAYS ###

static bool k_fill(vial_kctx_t* k, double v, vial_val_t* out) {
    const vial_type_t* t = k->type;
    uint32_t dims[VIAL_MAX_DIMS];
    for (uint8_t d = 0; d < t->ndim; ++d) {
        if (t->dims[d].kind != VIAL_EXT_STATIC) return vial_rt_fail(k, "R0103", "a shape not known when compiling");
        dims[d] = t->dims[d].n;
    }
    vial_val_t o; if (!vial_k_make(k, t->dtype, t->ndim, dims, &o)) return false;
    const size_t n = vial_val_count(&o);
    for (size_t i = 0; i < n; ++i) {
        switch (t->dtype) {
        case VIAL_DT_BOOL:    ((uint8_t*)o.data)[i] = v != 0.0; break;
        case VIAL_DT_INT:     ((int32_t*)o.data)[i] = (int32_t)v; break;
        case VIAL_DT_FLOAT:   ((float*)o.data)[i] = (float)v; break;
        case VIAL_DT_FLOAT64: ((double*)o.data)[i] = v; break;
        }
    }
    *out = o;
    return true;
}

static bool k_like(vial_kctx_t* k, const vial_val_t* a, double v, vial_val_t* out) {
    if (!a) return vial_rt_fail(k, "R0103", "a missing argument");
    vial_val_t o; if (!vial_k_make(k, k->type->dtype, a->ndim, a->dims, &o)) return false;
    const size_t n = vial_val_count(&o);
    for (size_t i = 0; i < n; ++i) {
        switch (o.dtype) {
        case VIAL_DT_BOOL:    ((uint8_t*)o.data)[i] = v != 0.0; break;
        case VIAL_DT_INT:     ((int32_t*)o.data)[i] = (int32_t)v; break;
        case VIAL_DT_FLOAT:   ((float*)o.data)[i] = (float)v; break;
        case VIAL_DT_FLOAT64: ((double*)o.data)[i] = v; break;
        }
    }
    *out = o;
    return true;
}

static bool k_range(vial_kctx_t* k, const bargs_t* ba, vial_val_t* out) {
    int64_t lo = 0, hi = 0, step = 1;
    const vial_val_t* a0 = arg(ba, "", 0);
    const vial_val_t* a1 = arg(ba, "", 1);
    const vial_val_t* a2 = arg(ba, "", 2);
    if (!a0) return vial_rt_fail(k, "R0103", "range takes 1 to 3 ints");
    if (a1) { lo = (int64_t)num_at(a0, 0); hi = (int64_t)num_at(a1, 0); } else hi = (int64_t)num_at(a0, 0);
    if (a2) step = (int64_t)num_at(a2, 0);
    if (!step) return vial_rt_fail(k, "R0103", "range's step is 0");
    int64_t c = step > 0 ? (hi - lo + step - 1) / step : (lo - hi - step - 1) / (-step);
    if (c < 0) c = 0;
    if (!vial_k_room(k, (size_t)c)) return false;
    const uint32_t n = (uint32_t)c;
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_INT, 1, &n, &o)) return false;
    for (uint32_t i = 0; i < n; ++i) ((int32_t*)o.data)[i] = (int32_t)(lo + (int64_t)i * step);
    *out = o;
    return true;
}

static bool k_concatenate(vial_kctx_t* k, const vial_val_t* items, const vial_val_t* axis, vial_val_t* out) {
    if (!items || items->dtype != VIAL_DT_TUPLE) return vial_rt_fail(k, "R0103", "concatenate takes a tuple or list of arrays");
    const uint32_t n = items->num_elems;
    if (!n) return vial_rt_fail(k, "R0103", "concatenate of nothing");
    const vial_val_t* e = items->elems;
    const uint8_t nd = e[0].ndim;
    int64_t ax = axis && axis->dtype == VIAL_DT_INT ? *(const int32_t*)axis->data : 0;
    if (ax < 0) ax += nd;
    if (ax < 0 || ax >= nd) return vial_rt_fail(k, "R0103", "axis %lld of %u dimensions", (long long)ax, nd);
    uint32_t dims[VIAL_MAX_DIMS];
    memcpy(dims, e[0].dims, nd * sizeof(uint32_t));
    dims[ax] = 0;
    for (uint32_t i = 0; i < n; ++i) {
        if (e[i].ndim != nd) return vial_rt_fail(k, "R0103", "concatenate of different ranks");
        for (uint8_t d = 0; d < nd; ++d) if (d != ax && e[i].dims[d] != e[0].dims[d]) return vial_rt_fail(k, "R0103", "concatenate of different shapes");
        dims[ax] += e[i].dims[ax];
    }
    const uint8_t dt = k->type->dtype;
    if (dt == VIAL_DT_SET) {
        if (nd != 1) return vial_rt_fail(k, "R0103", "concatenate of sets joins one dimension");
        return vial_k_stack(&(vial_kctx_t){.in = k->in, .op = k->op, .type = &(vial_type_t){.dtype = VIAL_DT_SET}, .alloc = k->alloc, .map_index = -1},
                            e, n, out) && (out->ndim = 1, out->dims[0] = dims[0], true);
    }
    vial_val_t o; if (!vial_k_make(k, dt, nd, dims, &o)) return false;
    const size_t es = vial_dtype_size(dt);
    size_t outer = 1, inner = 1;
    for (int64_t d = 0; d < ax; ++d) outer *= dims[d];
    for (uint8_t d = (uint8_t)(ax + 1); d < nd; ++d) inner *= dims[d];
    size_t at = 0;
    for (uint32_t i = 0; i < n; ++i) {
        const vial_val_t c = vial_val_as(&e[i], dt, k->alloc);
        const size_t chunk = (size_t)e[i].dims[ax] * inner;
        for (size_t q = 0; q < outer; ++q) {
            memcpy((uint8_t*)o.data + (q * dims[ax] * inner + at) * es, (const uint8_t*)c.data + q * chunk * es, chunk * es);
            for (size_t r = 0; r < chunk; ++r) if (vial_val_is_missing(&c, q * chunk + r)) vial_val_set_missing(&o, q * dims[ax] * inner + at + r, k->alloc);
        }
        at += chunk;
    }
    *out = o;
    return true;
}

static bool k_missing_where(vial_kctx_t* k, const vial_val_t* c, const vial_val_t* x, vial_val_t* out) {
    if (!c || !x) return vial_rt_fail(k, "R0103", "a missing argument");
    vial_val_t o = *x;
    o.missing = NULL;
    const size_t n = vial_val_count(x);
    // c broadcasts against x
    for (size_t i = 0; i < n; ++i) {
        size_t rem = i, off = 0, s = 1;
        for (int d = x->ndim - 1; d >= 0; --d) {
            const uint32_t q = (uint32_t)(rem % x->dims[d]);
            rem /= x->dims[d];
            const int dd = d - (x->ndim - c->ndim);
            if (dd >= 0) { if (c->dims[dd] != 1) off += q * s; s *= c->dims[dd]; }
        }
        if (vial_val_is_missing(x, i) || (((const uint8_t*)c->data)[off] && !vial_val_is_missing(c, off))) vial_val_set_missing(&o, i, k->alloc);
    }
    *out = o;
    return true;
}

static bool k_valid(vial_kctx_t* k, const vial_val_t* x, vial_val_t* out) {
    if (!x) return vial_rt_fail(k, "R0103", "a missing argument");
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_BOOL, x->ndim, x->dims, &o)) return false;
    const size_t n = vial_val_count(x);
    for (size_t i = 0; i < n; ++i) ((uint8_t*)o.data)[i] = !vial_val_is_missing(x, i);
    *out = o;
    return true;
}

static bool k_filled(vial_kctx_t* k, const vial_val_t* x, const vial_val_t* v, vial_val_t* out) {
    if (!x || !v) return vial_rt_fail(k, "R0103", "a missing argument");
    const uint8_t dt = k->type->dtype;
    vial_val_t o = vial_val_as(x, dt, k->alloc);
    vial_val_t c; if (!vial_k_make(k, dt, o.ndim, o.dims, &c)) return false;
    const size_t n = vial_val_count(&o), es = vial_dtype_size(dt);
    memcpy(c.data, o.data, n * es);
    const vial_val_t fv = vial_val_as(v, dt, k->alloc);
    for (size_t i = 0; i < n; ++i) if (vial_val_is_missing(x, i)) memcpy((uint8_t*)c.data + i * es, fv.data, es);
    *out = c;
    return true;
}

static bool k_diff(vial_kctx_t* k, const vial_val_t* a, const vial_val_t* axis, vial_val_t* out) {
    if (!a) return vial_rt_fail(k, "R0103", "a missing argument");
    const uint8_t dt = k->type->dtype;
    const vial_val_t x = vial_val_as(a, dt, k->alloc);
    int64_t ax = axis && axis->dtype == VIAL_DT_INT ? *(const int32_t*)axis->data : 0;
    if (ax < 0) ax += x.ndim;
    if (ax < 0 || ax >= x.ndim) return vial_rt_fail(k, "R0103", "axis out of range");
    uint32_t dims[VIAL_MAX_DIMS];
    memcpy(dims, x.dims, x.ndim * sizeof(uint32_t));
    if (!dims[ax]) return vial_rt_fail(k, "R0103", "diff of an empty dimension");
    dims[ax] -= 1;
    vial_val_t o; if (!vial_k_make(k, dt, x.ndim, dims, &o)) return false;
    size_t outer = 1, inner = 1;
    for (int64_t d = 0; d < ax; ++d) outer *= x.dims[d];
    for (uint8_t d = (uint8_t)(ax + 1); d < x.ndim; ++d) inner *= x.dims[d];
    for (size_t q = 0; q < outer; ++q) {
        for (uint32_t i = 0; i < dims[ax]; ++i) {
            for (size_t r = 0; r < inner; ++r) {
                const size_t s0 = (q * x.dims[ax] + i) * inner + r, s1 = s0 + inner, d0 = (q * dims[ax] + i) * inner + r;
                switch (dt) {
                case VIAL_DT_FLOAT:   ((float*)o.data)[d0] = ((const float*)x.data)[s1] - ((const float*)x.data)[s0]; break;
                case VIAL_DT_FLOAT64: ((double*)o.data)[d0] = ((const double*)x.data)[s1] - ((const double*)x.data)[s0]; break;
                case VIAL_DT_INT:     ((int32_t*)o.data)[d0] = ((const int32_t*)x.data)[s1] - ((const int32_t*)x.data)[s0]; break;
                }
                if (vial_val_is_missing(&x, s0) || vial_val_is_missing(&x, s1)) vial_val_set_missing(&o, d0, k->alloc);
            }
        }
    }
    *out = o;
    return true;
}

// ### MORE GEOMETRY ###

// Whether point a is within r of point b, d <= r with d as distance() gives it, bit for bit. A search's own squared
// distance decides where it is clearly on one side, beyond slack of r; between, the exact distance does.
typedef struct within_t {
    uint8_t*            hit;
    float               r;
    float               lo2, hi2;       // (r - slack)^2, -1 when r <= slack; (r + slack)^2
    const float*        p;              // the points of `of`
    const float*        xyz;            // every particle's position
    const md_unitcell_t* cell;
} within_t;

static inline bool within_decide(const within_t* w, float d2, uint32_t i, uint32_t j) {
    if (d2 <= w->lo2) return true;
    if (d2 > w->hi2) return false;
    return pair_entry(w->p + 3 * (size_t)i, w->xyz + 3 * (size_t)j, w->cell) <= w->r;   // near r, or not a number
}

static void within_cb(const uint32_t* i_idx, const uint32_t* j_idx, const float* d2, size_t n, void* user) {
    within_t* w = user;
    for (size_t q = 0; q < n; ++q) {
        const uint32_t j = j_idx[q];
        if (!w->hit[j] && within_decide(w, d2[q], i_idx[q], j)) w->hit[j] = 1;
    }
}

// Every particle against every point of `of`, eight particles at a time, skipping the groups of eight points of `of`
// too far from them: unless `of` is large, cheaper than building a search over every particle each frame. False
// where it does not apply: a triclinic cell.
static bool within_brute(vial_kctx_t* k, const within_t* w, size_t na, size_t np) {
    const md_unitcell_t* cell = w->cell;
    if (md_unitcell_flags(cell) & MD_UNITCELL_TRICLINIC) return false;
    float ext[3];
    md_unitcell_diag_extract_float(ext, cell);
    const bool ortho = (md_unitcell_flags(cell) & MD_UNITCELL_ORTHO) != 0;
    md_256 vext[3], vinv[3];
    for (int c = 0; c < 3; ++c) {
        const float e = ortho ? ext[c] : 0.0f;
        vext[c] = md_mm256_set1_ps(e);
        vinv[c] = md_mm256_set1_ps(e != 0.0f ? 1.0f / e : 0.0f);
    }
    // The particles, one array per axis, padded to eight with the last
    const size_t cap = (na + 7) & ~(size_t)7;
    float* ax = md_alloc(k->alloc, 3 * cap * sizeof(float));
    float* ay = ax + cap;
    float* az = ay + cap;
    for (size_t j = 0; j < cap; ++j) {
        const float* q = w->xyz + 3 * MIN(j, na - 1);
        ax[j] = q[0]; ay[j] = q[1]; az[j] = q[2];
    }
    // `of` in groups of eight points in order, each as one image around its first point: a centre and a radius. A
    // particle further from the centre than r, the radius and the rounding is within r of none of them.
    const size_t nc = (np + 7) / 8;
    float* cen = md_alloc(k->alloc, 4 * nc * sizeof(float));
    for (size_t c = 0; c < nc; ++c) {
        const float* p0 = w->p + 24 * c;
        const size_t cnt = MIN((size_t)8, np - 8 * c);
        float q[8][3], lo[3], hi[3];
        for (size_t i = 0; i < cnt; ++i) {
            for (int d = 0; d < 3; ++d) {
                float v = p0[3 * i + d] - p0[d];
                if (ortho && ext[d] != 0.0f) v -= ext[d] * nearbyintf(v / ext[d]);
                q[i][d] = p0[d] + v;
                lo[d] = i ? MIN(lo[d], q[i][d]) : q[i][d];
                hi[d] = i ? MAX(hi[d], q[i][d]) : q[i][d];
            }
        }
        float r2 = 0.0f;
        for (int d = 0; d < 3; ++d) cen[4 * c + d] = 0.5f * (lo[d] + hi[d]);
        for (size_t i = 0; i < cnt; ++i) {
            float e = 0.0f;
            for (int d = 0; d < 3; ++d) e += (q[i][d] - cen[4 * c + d]) * (q[i][d] - cen[4 * c + d]);
            r2 = MAX(r2, e);
        }
        cen[4 * c + 3] = sqrtf(r2);
    }
    const float margin = sqrtf(w->hi2) - w->r;   // the slack
    const md_256 lo2 = md_mm256_set1_ps(w->lo2), hi2 = md_mm256_set1_ps(w->hi2);
    for (size_t j0 = 0; j0 < na; j0 += 8) {
        const int valid = na - j0 >= 8 ? 0xFF : (1 << (na - j0)) - 1;
        const md_256 x = md_mm256_loadu_ps(ax + j0), y = md_mm256_loadu_ps(ay + j0), z = md_mm256_loadu_ps(az + j0);
        int hit = 0;
        for (size_t c = 0; c < nc && hit != valid; ++c) {
            md_256 dx = md_mm256_sub_ps(x, md_mm256_set1_ps(cen[4 * c + 0]));
            md_256 dy = md_mm256_sub_ps(y, md_mm256_set1_ps(cen[4 * c + 1]));
            md_256 dz = md_mm256_sub_ps(z, md_mm256_set1_ps(cen[4 * c + 2]));
            dx = md_mm256_sub_ps(dx, md_mm256_mul_ps(vext[0], md_mm256_round_ps(md_mm256_mul_ps(dx, vinv[0]))));
            dy = md_mm256_sub_ps(dy, md_mm256_mul_ps(vext[1], md_mm256_round_ps(md_mm256_mul_ps(dy, vinv[1]))));
            dz = md_mm256_sub_ps(dz, md_mm256_mul_ps(vext[2], md_mm256_round_ps(md_mm256_mul_ps(dz, vinv[2]))));
            const md_256 dc2 = md_mm256_add_ps(md_mm256_add_ps(md_mm256_mul_ps(dx, dx), md_mm256_mul_ps(dy, dy)), md_mm256_mul_ps(dz, dz));
            const float reach = w->r + cen[4 * c + 3] + 2.0f * margin;
            // Not a number passes, for the exact distance to decide
            const int near = (md_mm256_movemask_ps(md_mm256_cmpgt_ps(dc2, md_mm256_set1_ps(reach * reach))) ^ 0xFF) & valid & ~hit;
            if (!near) continue;
            for (size_t i = 8 * c; i < MIN(np, 8 * c + 8); ++i) {
                const float* p = w->p + 3 * i;
                dx = md_mm256_sub_ps(x, md_mm256_set1_ps(p[0]));
                dy = md_mm256_sub_ps(y, md_mm256_set1_ps(p[1]));
                dz = md_mm256_sub_ps(z, md_mm256_set1_ps(p[2]));
                dx = md_mm256_sub_ps(dx, md_mm256_mul_ps(vext[0], md_mm256_round_ps(md_mm256_mul_ps(dx, vinv[0]))));
                dy = md_mm256_sub_ps(dy, md_mm256_mul_ps(vext[1], md_mm256_round_ps(md_mm256_mul_ps(dy, vinv[1]))));
                dz = md_mm256_sub_ps(dz, md_mm256_mul_ps(vext[2], md_mm256_round_ps(md_mm256_mul_ps(dz, vinv[2]))));
                const md_256 d2 = md_mm256_add_ps(md_mm256_add_ps(md_mm256_mul_ps(dx, dx), md_mm256_mul_ps(dy, dy)), md_mm256_mul_ps(dz, dz));
                const int in = md_mm256_movemask_ps(md_mm256_cmple_ps(d2, lo2));
                const int out = md_mm256_movemask_ps(md_mm256_cmpgt_ps(d2, hi2));
                hit |= in & near;
                // Near r, or not a number: the exact distance decides
                int unsure = ~(in | out) & near & ~hit;
                while (unsure) {
                    const int lane = ctz32((uint32_t)unsure);
                    unsure &= unsure - 1;
                    if (pair_entry(p, w->xyz + 3 * (j0 + lane), cell) <= w->r) hit |= 1 << lane;
                }
            }
        }
        while (hit) {
            const int lane = ctz32((uint32_t)hit);
            hit &= hit - 1;
            w->hit[j0 + lane] = 1;
        }
    }
    return true;
}

// The scope's particles within r of any point of `of`, d <= r under the minimum image, d as distance() gives it
static bool k_within(vial_kctx_t* k, const vial_val_t* rv, const vial_val_t* of, vial_val_t* out) {
    const float* xyz = frame_xyz(k);
    if (!xyz) return false;
    if (!rv || !of) return vial_rt_fail(k, "R0103", "within takes r and of=");
    const float r = (float)num_at(rv, 0);
    vial_val_t p;
    if (!points(k, of, &p)) return false;
    const size_t np = vial_val_count(&p) / 3;
    const size_t na = k->in->frame->num_atoms;
    uint8_t* hit = md_alloc(k->alloc, na ? na : 1);
    memset(hit, 0, na);
    if (np && na && r >= 0.0f) {
        const md_unitcell_t cell = frame_cell(k);
        // How far a search's squared distance may be from the exact one's square (see pair_extreme_fast)
        float scale = 0.0f;
        finite_scale(xyz, 3 * na, &scale);
        finite_scale(p.data, 3 * np, &scale);
        float ext[3];
        md_unitcell_diag_extract_float(ext, &cell);
        const float slack = (scale + ext[0] + ext[1] + ext[2]) * 1.0e-5f + 1.0e-4f;
        within_t w = {.hit = hit, .r = r, .lo2 = r > slack ? (r - slack) * (r - slack) : -1.0f, .hi2 = (r + slack) * (r + slack),
                      .p = p.data, .xyz = xyz, .cell = &cell};
        // A few points of `of`: every pair, which costs less than building a search over every particle
        if (k->in->plain || ((np + 7) / 8) * ((na + 7) / 8) > VIAL_WITHIN_BRUTE || !within_brute(k, &w, na, np)) {
            md_coord_stream_t all = md_coord_stream_from_aos(xyz, 3 * sizeof(float), NULL, na);
            // The search may leave out a pair at exactly r; it is asked past r by the slack, and the exact distance
            // decides there. The structure is built for that same cutoff: built for r, a query past it reaches two
            // cells around, not one.
            const double reach = (double)r + slack;
            md_spatial_acc_t acc = {.alloc = k->alloc};
            k_spatial_init(k, &acc, &all, reach, &cell);
            md_coord_stream_t ext_pts = md_coord_stream_from_aos((const float*)p.data, 3 * sizeof(float), NULL, np);
            md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &ext_pts, reach, within_cb, &w, 0);
            md_spatial_acc_free(&acc);
        }
    }
    uint32_t n = 0;
    for (size_t i = 0; i < na; ++i) n += hit[i];
    uint32_t* idx = md_alloc(k->alloc, (n ? n : 1) * sizeof(uint32_t));
    uint32_t w = 0;
    for (size_t i = 0; i < na; ++i) if (hit[i]) idx[w++] = (uint32_t)i;
    vial_sets_t* s = vial_sets_one(idx, n, k->alloc);
    s->level = VIAL_LEVEL_PARTICLE;
    *out = vial_val_sets(s);
    return true;
}

// p superposed onto ref, by the rows of p whose particles are in `on` (all rows without it), with equal or mass weights
static bool k_fit(vial_kctx_t* k, const vial_val_t* pv, const vial_val_t* rv, const vial_val_t* on, const vial_val_t* weights, vial_val_t* out) {
    if (!pv || !rv) return vial_rt_fail(k, "R0103", "a missing argument");
    vial_val_t p, r;
    if (!points(k, pv, &p) || !points(k, rv, &r)) return false;
    if (p.ndim != 2 || r.ndim != 2) return vial_rt_fail(k, "R0103", "fit takes one group of points each");
    const uint32_t n = p.dims[0];
    if (r.dims[0] != n) return vial_rt_fail(k, "R0103", "fit pairs %u points with %u", n, r.dims[0]);
    const bool by_mass = str_is(weights, "mass");
    const md_system_t* sys = k->in->msys;
    // The rows that fit, and their weights
    uint8_t* use = md_alloc(k->alloc, n ? n : 1);
    double* w = md_alloc(k->alloc, (n ? n : 1) * sizeof(double));
    for (uint32_t i = 0; i < n; ++i) {
        const uint32_t particle = p.origin ? vial_sets_items(p.origin, 0)[i] : 0;
        use[i] = 1;
        if (on && on->dtype == VIAL_DT_SET) {
            if (!p.origin) return vial_rt_fail(k, "R0103", "fit on= needs points whose particles are known");
            use[i] = 0;
            const uint32_t m = vial_sets_size(on->sets, 0);
            const uint32_t* oi = vial_sets_items(on->sets, 0);
            for (uint32_t q = 0; q < m && !use[i]; ++q) use[i] = oi[q] == particle;
        }
        w[i] = 1.0;
        if (by_mass) {
            if (!p.origin || !sys) return vial_rt_fail(k, "R0103", "fit by mass needs points whose particles are known");
            w[i] = k_mass(k_masses(k), sys, particle);
        }
    }
    double cp[3] = {0}, cr[3] = {0}, ws = 0.0;
    for (uint32_t i = 0; i < n; ++i) {
        if (!use[i]) continue;
        for (int c = 0; c < 3; ++c) { cp[c] += w[i] * ((const float*)p.data)[3 * i + c]; cr[c] += w[i] * ((const float*)r.data)[3 * i + c]; }
        ws += w[i];
    }
    if (ws <= 0.0) return vial_rt_fail(k, "R0103", "fit on no particles");
    for (int c = 0; c < 3; ++c) { cp[c] /= ws; cr[c] /= ws; }
    vec4_t* rel_p = md_alloc(k->alloc, (n ? n : 1) * sizeof(vec4_t));
    vec4_t* rel_r = md_alloc(k->alloc, (n ? n : 1) * sizeof(vec4_t));
    uint32_t m = 0;
    for (uint32_t i = 0; i < n; ++i) {
        if (!use[i]) continue;
        const float* a = (const float*)p.data + 3 * i;
        const float* b = (const float*)r.data + 3 * i;
        rel_p[m] = vec4_set((float)(a[0] - cp[0]), (float)(a[1] - cp[1]), (float)(a[2] - cp[2]), (float)w[i]);
        rel_r[m] = vec4_set((float)(b[0] - cr[0]), (float)(b[1] - cr[1]), (float)(b[2] - cr[2]), (float)w[i]);
        m++;
    }
    const mat3_t R = md_util_optimal_rotation_rel_vec4(rel_r, rel_p, m);
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, 2, p.dims, &o)) return false;
    for (uint32_t i = 0; i < n; ++i) {
        const float* a = (const float*)p.data + 3 * i;
        const vec3_t v = mat3_mul_vec3(R, (vec3_t){{(float)(a[0] - cp[0]), (float)(a[1] - cp[1]), (float)(a[2] - cp[2])}});
        float* d = (float*)o.data + 3 * i;
        d[0] = (float)(v.x + cr[0]); d[1] = (float)(v.y + cr[1]); d[2] = (float)(v.z + cr[2]);
    }
    o.origin = p.origin;
    *out = o;
    return true;
}

// ### MORE ARRAYS ###

static bool k_component(vial_kctx_t* k, const vial_val_t* a, int c, vial_val_t* out) {
    if (!a) return vial_rt_fail(k, "R0103", "a missing argument");
    vial_val_t p;
    if (!points(k, a, &p)) return false;
    const size_t g = vial_val_count(&p) / 3;
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, p.ndim - 1, p.dims, &o)) return false;
    for (size_t i = 0; i < g; ++i) {
        ((float*)o.data)[i] = ((const float*)p.data)[3 * i + c];
        if (vial_val_is_missing(&p, 3 * i + c)) vial_val_set_missing(&o, i, k->alloc);
    }
    o.origin = p.origin;
    *out = o;
    return true;
}

static bool k_transpose(vial_kctx_t* k, const vial_val_t* a, vial_val_t* out) {
    if (!a) return vial_rt_fail(k, "R0103", "a missing argument");
    if (a->ndim < 2) { *out = *a; return true; }
    uint32_t dims[VIAL_MAX_DIMS];
    for (uint8_t d = 0; d < a->ndim; ++d) dims[d] = a->dims[a->ndim - 1 - d];
    vial_val_t o; if (!vial_k_make(k, a->dtype, a->ndim, dims, &o)) return false;
    const size_t n = vial_val_count(a), es = vial_dtype_size(a->dtype);
    size_t st[VIAL_MAX_DIMS];
    {
        size_t s = 1;
        for (int d = a->ndim - 1; d >= 0; --d) { st[d] = s; s *= a->dims[d]; }
    }
    for (size_t i = 0; i < n; ++i) {
        // Output index i over reversed dims: its source offset
        size_t rem = i, off = 0;
        for (int d = a->ndim - 1; d >= 0; --d) {
            const uint32_t x = (uint32_t)(rem % dims[d]);
            rem /= dims[d];
            off += x * st[a->ndim - 1 - d];
        }
        memcpy((uint8_t*)o.data + i * es, (const uint8_t*)a->data + off * es, es);
        if (vial_val_is_missing(a, off)) vial_val_set_missing(&o, i, k->alloc);
    }
    *out = o;
    return true;
}

static bool k_reshape(vial_kctx_t* k, const vial_val_t* a, vial_val_t* out) {
    if (!a) return vial_rt_fail(k, "R0103", "a missing argument");
    // The shape is the result type's; one extent may be known only now, as the rest's remainder
    const vial_type_t* t = k->type;
    vial_val_t o = *a;
    o.ndim = t->ndim;
    size_t known = 1;
    int open = -1;
    for (uint8_t d = 0; d < t->ndim; ++d) {
        if (t->dims[d].kind == VIAL_EXT_STATIC) { o.dims[d] = t->dims[d].n; known *= t->dims[d].n; }
        else if (open < 0) open = d;
        else return vial_rt_fail(k, "R0103", "a shape with two extents not known");
    }
    const size_t n = vial_val_count(a);
    if (open >= 0) { if (!known || n % known) return vial_rt_fail(k, "R0103", "%zu values do not reshape", n); o.dims[open] = (uint32_t)(n / known); }
    else if (known != n) return vial_rt_fail(k, "R0103", "%zu values do not reshape to %zu", n, known);
    o.origin = NULL;
    *out = o;
    return true;
}

static bool k_linspace(vial_kctx_t* k, const vial_val_t* lo, const vial_val_t* hi, const vial_val_t* nv, vial_val_t* out) {
    if (!lo || !hi || !nv) return vial_rt_fail(k, "R0103", "a missing argument");
    const double a = num_at(lo, 0), b = num_at(hi, 0);
    const double nd = num_at(nv, 0);
    if (!(nd >= 0.0) || !vial_k_room(k, nd > 1e18 ? SIZE_MAX : (size_t)nd)) return nd >= 0.0 ? false : vial_rt_fail(k, "R0103", "linspace of a negative count");
    const int64_t n = (int64_t)nd;
    const uint32_t cnt = (uint32_t)n;
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, 1, &cnt, &o)) return false;
    // NumPy's: start + i * step, and the last exactly hi
    const double step = cnt > 1 ? (b - a) / (double)(cnt - 1) : 0.0;
    for (uint32_t i = 0; i < cnt; ++i) ((float*)o.data)[i] = (float)(a + (double)i * step);
    if (cnt > 1) ((float*)o.data)[cnt - 1] = (float)b;
    *out = o;
    return true;
}

static bool k_arange(vial_kctx_t* k, const vial_val_t* lo, const vial_val_t* hi, const vial_val_t* sv, vial_val_t* out) {
    if (!lo || !hi) return vial_rt_fail(k, "R0103", "a missing argument");
    const double a = num_at(lo, 0), b = num_at(hi, 0), s = sv ? num_at(sv, 0) : 1.0;
    if (s == 0.0) return vial_rt_fail(k, "R0103", "arange's step is 0");
    double c = ceil((b - a) / s);
    if (!(c >= 0)) c = 0;
    if (!vial_k_room(k, c > 1e18 ? SIZE_MAX : (size_t)c)) return false;
    const uint32_t cnt = (uint32_t)c;
    vial_val_t o; if (!vial_k_make(k, k->type->dtype, 1, &cnt, &o)) return false;
    for (uint32_t i = 0; i < cnt; ++i) {
        const double v = a + (double)i * s;
        if (o.dtype == VIAL_DT_INT) ((int32_t*)o.data)[i] = (int32_t)v; else ((float*)o.data)[i] = (float)v;
    }
    *out = o;
    return true;
}

static bool k_eye(vial_kctx_t* k, const vial_val_t* nv, vial_val_t* out) {
    if (!nv) return vial_rt_fail(k, "R0103", "a missing argument");
    const double nd = num_at(nv, 0);
    if (!(nd >= 0.0) || nd > 65536.0) return vial_rt_fail(k, "R0105", "eye of %g", nd);
    const uint32_t n = (uint32_t)nd;
    vial_val_t o; if (!vial_k_make_zero(k, VIAL_DT_FLOAT, 2, (uint32_t[]){n, n}, &o)) return false;
    for (uint32_t i = 0; i < n; ++i) ((float*)o.data)[i * n + i] = 1.0f;
    *out = o;
    return true;
}

// Each particle of the groups gets its group's value; the result runs over the groups' particles, merged and ascending
static bool k_spread(vial_kctx_t* k, const vial_val_t* values, const vial_val_t* groups, const vial_val_t* combine, vial_val_t* out) {
    if (!values) return vial_rt_fail(k, "R0103", "a missing argument");
    if (!groups || groups->dtype != VIAL_DT_SET || groups->ndim != 1) return vial_rt_fail(k, "R0103", "spread takes groups");
    const uint32_t g = groups->dims[0];
    if (!values->ndim || values->dims[0] != g) return vial_rt_fail(k, "R0103", "spread gives each of %u groups one value", g);
    const vial_sets_t* merged = vial_sets_merge(groups->sets, k->alloc);
    const uint32_t n = vial_sets_size(merged, 0);
    const uint32_t* mi = vial_sets_items(merged, 0);
    const uint8_t dt = k->type->dtype;
    const vial_val_t v = vial_val_as(values, dt == VIAL_DT_BOOL ? VIAL_DT_BOOL : VIAL_DT_FLOAT64, k->alloc);
    size_t inner = 1;
    for (uint8_t d = 1; d < values->ndim; ++d) inner *= values->dims[d];
    uint32_t dims[VIAL_MAX_DIMS];
    memcpy(dims, values->dims, values->ndim * sizeof(uint32_t));
    dims[0] = n;
    double* acc = md_alloc(k->alloc, ((size_t)n * inner ? (size_t)n * inner : 1) * sizeof(double));
    uint32_t* cnt = md_alloc(k->alloc, (n ? n : 1) * sizeof(uint32_t));
    memset(acc, 0, (size_t)n * inner * sizeof(double));
    memset(cnt, 0, n * sizeof(uint32_t));
    const bool mean = str_is(combine, "mean"), mx = str_is(combine, "max");
    for (uint32_t q = 0; q < g; ++q) {
        const uint32_t sz = vial_sets_size(groups->sets, q);
        const uint32_t* gi = vial_sets_items(groups->sets, q);
        for (uint32_t a = 0; a < sz; ++a) {
            // The particle's place among the merged ones
            uint32_t lo = 0, hi = n;
            while (lo < hi) { const uint32_t mid = (lo + hi) / 2; if (mi[mid] < gi[a]) lo = mid + 1; else hi = mid; }
            for (size_t r = 0; r < inner; ++r) {
                const double x = v.dtype == VIAL_DT_BOOL ? ((const uint8_t*)v.data)[q * inner + r] : ((const double*)v.data)[q * inner + r];
                double* dst = &acc[(size_t)lo * inner + r];
                if (!cnt[lo]) *dst = x;
                else if (mean) *dst += x;
                else if (mx) *dst = x > *dst ? x : *dst;
            }
            cnt[lo] += 1;
        }
    }
    vial_val_t o; if (!vial_k_make(k, dt, values->ndim, dims, &o)) return false;
    for (uint32_t p = 0; p < n; ++p) {
        for (size_t r = 0; r < inner; ++r) {
            double x = acc[(size_t)p * inner + r];
            if (mean && cnt[p]) x /= cnt[p];
            switch (dt) {
            case VIAL_DT_BOOL:    ((uint8_t*)o.data)[p * inner + r] = x != 0.0; break;
            case VIAL_DT_INT:     ((int32_t*)o.data)[p * inner + r] = (int32_t)x; break;
            case VIAL_DT_FLOAT:   ((float*)o.data)[p * inner + r] = (float)x; break;
            case VIAL_DT_FLOAT64: ((double*)o.data)[p * inner + r] = x; break;
            }
        }
    }
    o.origin = merged;
    *out = o;
    return true;
}

// ### MAPPING ###
// Built-ins whose kernels take one call's arguments are mapped here over the dimensions their kinds leave, as a def
// with kinds is: closest_pair, farthest_pair and fit.

static uint8_t core_rank(uint8_t kind, const vial_val_t* v) {
    switch (kind) {
    case VX_KIND_POINT:  return v->dtype == VIAL_DT_SET ? 0 : 1;
    case VX_KIND_POINTS: return v->dtype == VIAL_DT_SET ? 0 : (v->ndim >= 2 ? 2 : 1);
    case VX_KIND_GROUP:  return v->dtype == VIAL_DT_SET ? 0 : 1;
    case VX_KIND_SET: case VX_KIND_FLOAT: case VX_KIND_INT: case VX_KIND_BOOL: return 0;
    default:             return v->ndim;
    }
}

typedef bool (*core_fn)(vial_kctx_t* k, const vial_val_t* const* a, vial_val_t* out);

static bool distinct_of(const vial_val_t* d) { return d && d->dtype == VIAL_DT_BOOL && d->data && *(const uint8_t*)d->data; }
static bool core_pair_c(vial_kctx_t* k, const vial_val_t* const* a, vial_val_t* out) {
    return k_pair_x(k, a[0], a[1], true, distinct_of(a[2]), k->map_index, out);
}
static bool core_pair_f(vial_kctx_t* k, const vial_val_t* const* a, vial_val_t* out) {
    return k_pair_x(k, a[0], a[1], false, distinct_of(a[2]), k->map_index, out);
}
static bool core_fit(vial_kctx_t* k, const vial_val_t* const* a, vial_val_t* out) { return k_fit(k, a[0], a[1], a[2], a[3], out); }

static vial_val_t elem_of(vial_kctx_t* k, const vial_val_t* a, uint8_t lead, size_t i) {
    vial_val_t e = *a;
    e.ndim = a->ndim - lead;
    memmove(e.dims, a->dims + lead, e.ndim * sizeof(uint32_t));
    size_t inner = 1;
    for (uint8_t q = 0; q < e.ndim; ++q) inner *= e.dims[q];
    if (a->dtype == VIAL_DT_SET) {
        uint32_t pos = (uint32_t)i;
        vial_sets_t* s = vial_sets_pick(a->sets, &pos, 1, k->alloc);
        s->single = true;
        e.sets = s;
    } else if (a->data) {
        e.data = (uint8_t*)a->data + i * inner * vial_dtype_size(a->dtype);
        e.missing = NULL;
        for (size_t q = 0; q < inner; ++q) if (vial_val_is_missing(a, i * inner + q)) vial_val_set_missing(&e, q, k->alloc);
    }
    return e;
}

static bool mapped(vial_kctx_t* k, const vial_val_t* const* args, const uint8_t* kinds, uint32_t n, core_fn fn, vial_val_t* out) {
    uint8_t lead[8], lnd = 0;
    for (uint32_t i = 0; i < n; ++i) {
        lead[i] = args[i] ? (uint8_t)(args[i]->ndim - core_rank(kinds[i], args[i])) : 0;
        if (args[i] && args[i]->dtype == VIAL_DT_SET && args[i]->ndim == 0) lead[i] = 0;
        if (lead[i] > lnd) lnd = lead[i];
    }
    if (!lnd) return fn(k, args, out);
    uint32_t ldims[VIAL_MAX_DIMS];
    for (uint8_t d = 0; d < lnd; ++d) ldims[d] = 1;
    for (uint32_t i = 0; i < n; ++i) {
        for (uint8_t d = 0; d < lead[i]; ++d) {
            uint32_t* o = &ldims[lnd - lead[i] + d];
            const uint32_t e = args[i]->dims[d];
            if (*o == 1) *o = e;
            else if (e != 1 && e != *o) return vial_rt_fail(k, "R0103", "mapped extents %u and %u do not broadcast", *o, e);
        }
    }
    size_t count = 1;
    for (uint8_t d = 0; d < lnd; ++d) count *= ldims[d];
    vial_val_t* res = md_alloc(k->alloc, (count ? count : 1) * sizeof(vial_val_t));
    for (size_t c = 0; c < count; ++c) {
        uint32_t idx[VIAL_MAX_DIMS];
        size_t rem = c;
        for (int d = lnd - 1; d >= 0; --d) { idx[d] = (uint32_t)(rem % ldims[d]); rem /= ldims[d]; }
        vial_val_t el[8];
        const vial_val_t* ep[8];
        for (uint32_t i = 0; i < n; ++i) {
            if (!args[i] || !lead[i]) { ep[i] = args[i]; continue; }
            size_t lin = 0;
            for (uint8_t d = 0; d < lead[i]; ++d) {
                const uint32_t e = args[i]->dims[d];
                lin = lin * e + (e == 1 ? 0 : idx[lnd - lead[i] + d]);
            }
            el[i] = elem_of(k, args[i], lead[i], lin);
            ep[i] = &el[i];
        }
        k->map_index = lnd == 1 ? (int64_t)c : -1;
        const bool ok = fn(k, ep, &res[c]);
        k->map_index = -1;
        if (!ok) return false;
    }
    // Stacked: a tuple's elements each stacked, or the values themselves
    if (res[0].dtype == VIAL_DT_TUPLE) {
        const uint32_t ne = res[0].num_elems;
        vial_val_t* el = md_alloc(k->alloc, ne * sizeof(vial_val_t));
        for (uint32_t e = 0; e < ne; ++e) {
            const vial_val_t* f = &res[0].elems[e];
            uint32_t dims[VIAL_MAX_DIMS];
            memcpy(dims, ldims, lnd * sizeof(uint32_t));
            memcpy(dims + lnd, f->dims, f->ndim * sizeof(uint32_t));
            if (!vial_k_make(k, f->dtype, lnd + f->ndim, dims, &el[e])) return false;
            const size_t inner = vial_val_count(f), es = vial_dtype_size(f->dtype);
            for (size_t c = 0; c < count; ++c) {
                memcpy((uint8_t*)el[e].data + c * inner * es, res[c].elems[e].data, inner * es);
                for (size_t q = 0; q < inner; ++q) if (vial_val_is_missing(&res[c].elems[e], q)) vial_val_set_missing(&el[e], c * inner + q, k->alloc);
            }
        }
        *out = (vial_val_t){.dtype = VIAL_DT_TUPLE, .elems = el, .num_elems = ne};
        return true;
    }
    const vial_val_t* f = &res[0];
    uint32_t dims[VIAL_MAX_DIMS];
    memcpy(dims, ldims, lnd * sizeof(uint32_t));
    memcpy(dims + lnd, f->dims, f->ndim * sizeof(uint32_t));
    vial_val_t o; if (!vial_k_make(k, f->dtype, lnd + f->ndim, dims, &o)) return false;
    const size_t inner = vial_val_count(f), es = vial_dtype_size(f->dtype);
    for (size_t c = 0; c < count; ++c) {
        if (vial_val_count(&res[c]) != inner) return vial_rt_fail(k, "R0103", "mapped calls give values of different shapes");
        memcpy((uint8_t*)o.data + c * inner * es, res[c].data, inner * es);
        for (size_t q = 0; q < inner; ++q) if (vial_val_is_missing(&res[c], q)) vial_val_set_missing(&o, c * inner + q, k->alloc);
    }
    *out = o;
    return true;
}

// ### GROUPS: by= ###

// A reduction in groups along a's dimension over particles: one result per group, over the group's particles there
static bool k_reduce_by(vial_kctx_t* k, red_t red, const vial_val_t* a, const vial_val_t* by, vial_val_t* out) {
    if (!a) return vial_rt_fail(k, "R0103", "a missing argument");
    if (!by || by->dtype != VIAL_DT_SET || by->ndim != 1) return vial_rt_fail(k, "R0103", "by= takes an array of sets");
    // The particles a's first dimension runs over: from the values, or the type's description
    const vial_sets_t* over = a->origin;
    const vial_type_t* at = &k->in->ir->values[k->in->ir->args[k->op->arg_beg]].type;
    if (!over && at->ndim && at->dims[0].desc && k->in->ir->descs) over = k->in->ir->descs[at->dims[0].desc].sets;
    if (!over || !a->ndim || vial_sets_size(over, 0) != a->dims[0]) return vial_rt_fail(k, "R0103", "by= on values whose particles are not known");
    const uint32_t g = by->dims[0], n = a->dims[0];
    const uint32_t* oi = vial_sets_items(over, 0);
    size_t inner = 1;
    for (uint8_t d = 1; d < a->ndim; ++d) inner *= a->dims[d];
    // One value per group: every value of the group's rows, all their dimensions, as a reduction without axis= takes
    vial_val_t o; if (!vial_k_make(k, k->type->dtype, 1, &g, &o)) return false;
    vial_kctx_t k2 = *k;
    k2.op = NULL;
    vial_type_t t1 = *k->type;
    t1.ndim = 0;
    k2.type = &t1;
    const size_t es = vial_dtype_size(a->dtype), os = vial_dtype_size(o.dtype);
    for (uint32_t q = 0; q < g; ++q) {
        // The group's rows of a, in a's order
        const uint32_t sz = vial_sets_size(by->sets, q);
        const uint32_t* gi = vial_sets_items(by->sets, q);
        uint32_t* rows = md_alloc(k->alloc, (n ? n : 1) * sizeof(uint32_t));
        uint32_t nr = 0;
        for (uint32_t r = 0; r < n; ++r) {
            uint32_t lo = 0, hi = sz;
            while (lo < hi) { const uint32_t mid = (lo + hi) / 2; if (gi[mid] < oi[r]) lo = mid + 1; else hi = mid; }
            if (lo < sz && gi[lo] == oi[r]) rows[nr++] = r;
        }
        uint32_t sdims[VIAL_MAX_DIMS];
        memcpy(sdims, a->dims, a->ndim * sizeof(uint32_t));
        sdims[0] = nr;
        vial_val_t sub; if (!vial_k_make(k, a->dtype, a->ndim, sdims, &sub)) return false;
        for (uint32_t r = 0; r < nr; ++r) {
            memcpy((uint8_t*)sub.data + r * inner * es, (const uint8_t*)a->data + rows[r] * inner * es, inner * es);
            for (size_t z = 0; z < inner; ++z) if (vial_val_is_missing(a, rows[r] * inner + z)) vial_val_set_missing(&sub, r * inner + z, k->alloc);
        }
        vial_val_t one;
        if (!k_reduce(&k2, red, &sub, NULL, &one)) return false;
        memcpy((uint8_t*)o.data + q * os, one.data, os);
        if (vial_val_is_missing(&one, 0)) vial_val_set_missing(&o, q, k->alloc);
    }
    *out = o;
    return true;
}

// ### MATRICES ###

// NumPy's matmul: a vector is a row on the left and a column on the right, and leading dimensions broadcast
static bool k_matmul(vial_kctx_t* k, const vial_val_t* a, const vial_val_t* b, vial_val_t* out) {
    if (!a || !b) return vial_rt_fail(k, "R0103", "a missing argument");
    const vial_val_t x = vial_val_as(a, VIAL_DT_FLOAT, k->alloc), y = vial_val_as(b, VIAL_DT_FLOAT, k->alloc);
    if (!x.ndim || !y.ndim) return vial_rt_fail(k, "R0103", "@ takes arrays");
    const bool xv = x.ndim == 1, yv = y.ndim == 1;
    const uint32_t n = xv ? 1 : x.dims[x.ndim - 2], kk = x.dims[x.ndim - 1];
    const uint32_t kk2 = yv ? y.dims[0] : y.dims[y.ndim - 2], m = yv ? 1 : y.dims[y.ndim - 1];
    if (kk != kk2) return vial_rt_fail(k, "R0103", "@ of %u columns and %u rows", kk, kk2);
    const uint8_t bx = xv ? 0 : x.ndim - 2, by = yv ? 0 : y.ndim - 2;
    const uint8_t nb = bx > by ? bx : by;
    uint32_t bd[VIAL_MAX_DIMS];
    for (uint8_t d = 0; d < nb; ++d) bd[d] = 1;
    for (int w = 0; w < 2; ++w) {
        const vial_val_t* v = w ? &y : &x;
        const uint8_t nv = w ? by : bx;
        for (uint8_t d = 0; d < nv; ++d) {
            uint32_t* o = &bd[nb - nv + d];
            if (*o == 1) *o = v->dims[d];
            else if (v->dims[d] != 1 && v->dims[d] != *o) return vial_rt_fail(k, "R0103", "@'s leading dimensions do not broadcast");
        }
    }
    uint32_t od[VIAL_MAX_DIMS];
    uint8_t ond = 0;
    for (uint8_t d = 0; d < nb; ++d) od[ond++] = bd[d];
    if (!xv) od[ond++] = n;
    if (!yv) od[ond++] = m;
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, ond, od, &o)) return false;
    size_t nbatch = 1;
    for (uint8_t d = 0; d < nb; ++d) nbatch *= bd[d];
    for (size_t q = 0; q < nbatch; ++q) {
        size_t rem = q, ox = 0, oy = 0, sx = 1, sy = 1;
        for (int d = nb - 1; d >= 0; --d) {
            const uint32_t i = (uint32_t)(rem % bd[d]);
            rem /= bd[d];
            const int dx = d - (nb - bx), dy = d - (nb - by);
            if (dx >= 0) { if (x.dims[dx] != 1) ox += i * sx; sx *= x.dims[dx]; }
            if (dy >= 0) { if (y.dims[dy] != 1) oy += i * sy; sy *= y.dims[dy]; }
        }
        const float* A = (const float*)x.data + ox * n * kk;
        const float* B = (const float*)y.data + oy * kk * m;
        float* C = (float*)o.data + q * n * m;
        for (uint32_t i = 0; i < n; ++i) {
            for (uint32_t j = 0; j < m; ++j) {
                float s = 0.0f;
                for (uint32_t z = 0; z < kk; ++z) s += A[i * kk + z] * B[z * m + j];
                C[i * m + j] = s;
            }
        }
    }
    *out = o;
    return true;
}

bool vial_k_matmul(vial_kctx_t* k, const vial_val_t* a, const vial_val_t* b, vial_val_t* out) { return k_matmul(k, a, b, out); }

// Square matrices over the last two dimensions, each in doubles: n, and how many
static bool square(vial_kctx_t* k, const vial_val_t* m, uint32_t* n, size_t* count) {
    if (m->ndim < 2 || m->dims[m->ndim - 1] != m->dims[m->ndim - 2]) return vial_rt_fail(k, "R0103", "a square matrix is [n][n]");
    *n = m->dims[m->ndim - 1];
    *count = *n ? vial_val_count(m) / ((size_t)*n * *n) : 0;
    return true;
}

// Symmetric eigendecomposition: Householder tridiagonalization and the implicit QL algorithm (tred2 / tql2 as in
// EISPACK and JAMA; the same as md_gto_iso.c's), in doubles. In: V = S (row major). Out: V's columns the eigenvectors,
// val the eigenvalues; e is scratch of n. About 9 n^3 operations against cyclic Jacobi's 3 n^3 per sweep, with
// several sweeps, and it reads rows rather than columns: a 433 x 433 matrix in milliseconds rather than seconds.
static void sym_eigen(double* V, double* val, int n, double* e) {
#define EIG_V(i, j) V[(size_t)(i) * n + (j)]
    double* d = val;
    for (int j = 0; j < n; ++j) e[j] = 0.0;
    for (int j = 0; j < n; ++j) d[j] = EIG_V(n - 1, j);
    for (int i = n - 1; i > 0; --i) {
        double scale = 0.0, h = 0.0;
        for (int k = 0; k < i; ++k) scale += fabs(d[k]);
        if (scale == 0.0) {
            e[i] = d[i - 1];
            for (int j = 0; j < i; ++j) { d[j] = EIG_V(i - 1, j); EIG_V(i, j) = 0.0; EIG_V(j, i) = 0.0; }
        } else {
            for (int k = 0; k < i; ++k) { d[k] /= scale; h += d[k] * d[k]; }
            double f = d[i - 1];
            double g = sqrt(h);
            if (f > 0) g = -g;
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;
            for (int j = 0; j < i; ++j) e[j] = 0.0;
            for (int j = 0; j < i; ++j) {
                f = d[j];
                EIG_V(j, i) = f;
                g = e[j] + EIG_V(j, j) * f;
                for (int k = j + 1; k <= i - 1; ++k) { g += EIG_V(k, j) * d[k]; e[k] += EIG_V(k, j) * f; }
                e[j] = g;
            }
            f = 0.0;
            for (int j = 0; j < i; ++j) { e[j] /= h; f += e[j] * d[j]; }
            const double hh = f / (h + h);
            for (int j = 0; j < i; ++j) e[j] -= hh * d[j];
            for (int j = 0; j < i; ++j) {
                f = d[j];
                g = e[j];
                for (int k = j; k <= i - 1; ++k) EIG_V(k, j) -= (f * e[k] + g * d[k]);
                d[j] = EIG_V(i - 1, j);
                EIG_V(i, j) = 0.0;
            }
        }
        d[i] = h;
    }
    for (int i = 0; i < n - 1; ++i) {
        EIG_V(n - 1, i) = EIG_V(i, i);
        EIG_V(i, i) = 1.0;
        const double h = d[i + 1];
        if (h != 0.0) {
            for (int k = 0; k <= i; ++k) d[k] = EIG_V(k, i + 1) / h;
            for (int j = 0; j <= i; ++j) {
                double g = 0.0;
                for (int k = 0; k <= i; ++k) g += EIG_V(k, i + 1) * EIG_V(k, j);
                for (int k = 0; k <= i; ++k) EIG_V(k, j) -= g * d[k];
            }
        }
        for (int k = 0; k <= i; ++k) EIG_V(k, i + 1) = 0.0;
    }
    for (int j = 0; j < n; ++j) { d[j] = EIG_V(n - 1, j); EIG_V(n - 1, j) = 0.0; }
    EIG_V(n - 1, n - 1) = 1.0;
    e[0] = 0.0;

    for (int i = 1; i < n; ++i) e[i - 1] = e[i];
    e[n - 1] = 0.0;
    double f = 0.0, tst1 = 0.0;
    const double eps = DBL_EPSILON;
    for (int l = 0; l < n; ++l) {
        tst1 = fmax(tst1, fabs(d[l]) + fabs(e[l]));
        int m = l;
        while (m < n) {
            if (fabs(e[m]) <= eps * tst1) break;
            m++;
        }
        if (m > l) {
            int iter = 0;
            do {
                if (++iter > 60) break;
                double g = d[l];
                double p = (d[l + 1] - g) / (2.0 * e[l]);
                double r = hypot(p, 1.0);
                if (p < 0) r = -r;
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                const double dl1 = d[l + 1];
                double h = g - d[l];
                for (int i = l + 2; i < n; ++i) d[i] -= h;
                f += h;
                p = d[m];
                double c = 1.0, c2 = c, c3 = c;
                const double el1 = e[l + 1];
                double s = 0.0, s2 = 0.0;
                for (int i = m - 1; i >= l; --i) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypot(p, e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    for (int k = 0; k < n; ++k) {
                        h = EIG_V(k, i + 1);
                        EIG_V(k, i + 1) = s * EIG_V(k, i) + c * h;
                        EIG_V(k, i) = c * EIG_V(k, i) - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
            } while (fabs(e[l]) > eps * tst1);
        }
        d[l] += f;
        e[l] = 0.0;
    }
#undef EIG_V
}

static bool k_eigh(vial_kctx_t* k, const vial_val_t* m, vial_val_t* out) {
    if (!m) return vial_rt_fail(k, "R0103", "a missing argument");
    uint32_t n;
    size_t cnt;
    if (!square(k, m, &n, &cnt)) return false;
    vial_val_t w; if (!vial_k_make(k, VIAL_DT_FLOAT, m->ndim - 1, m->dims, &w)) return false;
    vial_val_t v; if (!vial_k_make(k, VIAL_DT_FLOAT, m->ndim, m->dims, &v)) return false;
    double* A = md_alloc(k->alloc, ((size_t)n * n ? (size_t)n * n : 1) * sizeof(double));
    double* V = md_alloc(k->alloc, ((size_t)n * n ? (size_t)n * n : 1) * sizeof(double));
    uint32_t* order = md_alloc(k->alloc, (n ? n : 1) * sizeof(uint32_t));
    double* val = md_alloc(k->alloc, (n ? n : 1) * sizeof(double));
    double* scratch = md_alloc(k->alloc, (n ? n : 1) * sizeof(double));
    for (size_t q = 0; q < cnt; ++q) {
        // The lower triangle, as NumPy's eigh reads it
        for (uint32_t i = 0; i < n; ++i) for (uint32_t j = 0; j < n; ++j) {
            const uint32_t a = i >= j ? i : j, b = i >= j ? j : i;
            A[i * n + j] = num_at(m, q * n * n + (size_t)a * n + b);
        }
        // Eigenvalues ascending, eigenvectors in the columns, each with its largest component positive so the sign
        // is the same on every machine
        memcpy(V, A, (size_t)n * n * sizeof(double));
        if (n) sym_eigen(V, val, (int)n, scratch);
        for (uint32_t i = 0; i < n; ++i) order[i] = i;
        for (uint32_t i = 1; i < n; ++i) {
            const uint32_t x = order[i];
            uint32_t j = i;
            while (j > 0 && val[order[j - 1]] > val[x]) { order[j] = order[j - 1]; --j; }
            order[j] = x;
        }
        for (uint32_t c = 0; c < n; ++c) {
            const uint32_t src = order[c];
            ((float*)w.data)[q * n + c] = (float)val[src];
            uint32_t big = 0;
            for (uint32_t r = 1; r < n; ++r) if (fabs(V[r * n + src]) > fabs(V[big * n + src])) big = r;
            const double sign = V[big * n + src] < 0 ? -1.0 : 1.0;
            for (uint32_t r = 0; r < n; ++r) ((float*)v.data)[q * n * n + (size_t)r * n + c] = (float)(sign * V[r * n + src]);
        }
    }
    vial_val_t* el = md_alloc(k->alloc, 2 * sizeof(vial_val_t));
    el[0] = w;
    el[1] = v;
    *out = (vial_val_t){.dtype = VIAL_DT_TUPLE, .elems = el, .num_elems = 2};
    return true;
}

// LU with partial pivoting, in doubles: the determinant, and solutions for rhs columns
static bool lu_solve(double* A, double* B, uint32_t n, uint32_t nrhs, double* det) {
    double d = 1.0;
    for (uint32_t c = 0; c < n; ++c) {
        uint32_t p = c;
        for (uint32_t r = c + 1; r < n; ++r) if (fabs(A[r * n + c]) > fabs(A[p * n + c])) p = r;
        if (A[p * n + c] == 0.0) { *det = 0.0; return false; }
        if (p != c) {
            for (uint32_t j = 0; j < n; ++j) { const double t = A[c * n + j]; A[c * n + j] = A[p * n + j]; A[p * n + j] = t; }
            for (uint32_t j = 0; j < nrhs; ++j) { const double t = B[c * nrhs + j]; B[c * nrhs + j] = B[p * nrhs + j]; B[p * nrhs + j] = t; }
            d = -d;
        }
        d *= A[c * n + c];
        for (uint32_t r = c + 1; r < n; ++r) {
            const double f = A[r * n + c] / A[c * n + c];
            for (uint32_t j = c; j < n; ++j) A[r * n + j] -= f * A[c * n + j];
            for (uint32_t j = 0; j < nrhs; ++j) B[r * nrhs + j] -= f * B[c * nrhs + j];
        }
    }
    for (int r = (int)n - 1; r >= 0; --r) {
        for (uint32_t j = 0; j < nrhs; ++j) {
            double s = B[r * nrhs + j];
            for (uint32_t c = r + 1; c < n; ++c) s -= A[r * n + c] * B[c * nrhs + j];
            B[r * nrhs + j] = s / A[r * n + r];
        }
    }
    *det = d;
    return true;
}

static bool k_linalg(vial_kctx_t* k, const char* name, const vial_val_t* m, const vial_val_t* b, vial_val_t* out) {
    if (!m) return vial_rt_fail(k, "R0103", "a missing argument");
    uint32_t n;
    size_t cnt;
    if (!square(k, m, &n, &cnt)) return false;
    const bool is_trace = !strcmp(name, "trace"), is_det = !strcmp(name, "det"), is_inv = !strcmp(name, "inv");
    vial_val_t o;
    uint32_t nrhs = n;
    if (is_trace || is_det) { if (!vial_k_make(k, VIAL_DT_FLOAT, m->ndim - 2, m->dims, &o)) return false; }
    else if (is_inv) { if (!vial_k_make(k, VIAL_DT_FLOAT, m->ndim, m->dims, &o)) return false; }
    else {
        // solve(m, b): b a vector [n] or a matrix [n][k]
        if (!b || !b->ndim || b->dims[0] != n || cnt != 1) return vial_rt_fail(k, "R0103", "solve takes [n][n] and [n] or [n][k]");
        nrhs = b->ndim == 1 ? 1 : b->dims[1];
        if (!vial_k_make(k, VIAL_DT_FLOAT, b->ndim, b->dims, &o)) return false;
    }
    double* A = md_alloc(k->alloc, ((size_t)n * n ? (size_t)n * n : 1) * sizeof(double));
    double* B = md_alloc(k->alloc, ((size_t)n * nrhs ? (size_t)n * nrhs : 1) * sizeof(double));
    for (size_t q = 0; q < cnt; ++q) {
        for (size_t i = 0; i < (size_t)n * n; ++i) A[i] = num_at(m, q * n * n + i);
        if (is_trace) {
            double t = 0.0;
            for (uint32_t i = 0; i < n; ++i) t += A[i * n + i];
            ((float*)o.data)[q] = (float)t;
            continue;
        }
        if (is_inv) for (uint32_t i = 0; i < n; ++i) for (uint32_t j = 0; j < n; ++j) B[i * n + j] = i == j;
        else if (!is_det) for (size_t i = 0; i < (size_t)n * nrhs; ++i) B[i] = num_at(b, i);
        double det;
        const bool ok = lu_solve(A, B, n, is_det ? 0 : nrhs, &det);
        if (is_det) { ((float*)o.data)[q] = (float)det; continue; }
        const size_t sz = (size_t)n * nrhs;
        for (size_t i = 0; i < sz; ++i) {
            ((float*)o.data)[q * sz + i] = (float)B[i];
            if (!ok) vial_val_set_missing(&o, q * sz + i, k->alloc);     // singular
        }
    }
    *out = o;
    return true;
}

// The covariance of the rows of x along axis: x[axis] are the samples, the rest the variables
static bool k_covariance(vial_kctx_t* k, const vial_val_t* x, const vial_val_t* axis, vial_val_t* out) {
    if (!x) return vial_rt_fail(k, "R0103", "a missing argument");
    if (!axis || axis->dtype == VIAL_DT_STR) return vial_rt_fail(k, "R0103", "covariance over frames is a trajectory value");
    if (x->ndim != 2) return vial_rt_fail(k, "R0103", "covariance takes samples by variables, [m][n]");
    int64_t ax = *(const int32_t*)axis->data;
    if (ax < 0) ax += 2;
    const uint32_t m = x->dims[ax], n = x->dims[1 - ax];
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, 2, (uint32_t[]){n, n}, &o)) return false;
    double* mean = md_alloc(k->alloc, (n ? n : 1) * sizeof(double));
    for (uint32_t j = 0; j < n; ++j) {
        double s = 0.0;
        for (uint32_t i = 0; i < m; ++i) s += num_at(x, ax == 0 ? (size_t)i * n + j : (size_t)j * m + i);
        mean[j] = m ? s / m : 0.0;
    }
    for (uint32_t a = 0; a < n; ++a) {
        for (uint32_t b = 0; b < n; ++b) {
            double s = 0.0;
            for (uint32_t i = 0; i < m; ++i) {
                const double xa = num_at(x, ax == 0 ? (size_t)i * n + a : (size_t)a * m + i) - mean[a];
                const double xb = num_at(x, ax == 0 ? (size_t)i * n + b : (size_t)b * m + i) - mean[b];
                s += xa * xb;
            }
            ((float*)o.data)[a * n + b] = (float)(m > 1 ? s / (m - 1) : 0.0);
        }
    }
    *out = o;
    return true;
}

// ### SELECTORS, WHEN EVALUATED ###
// What the check could not select when compiling: a positional set that changes per frame, as residues(within(...)).
// Keywords are the same filters the check builds.

static uint8_t level_of(const char* name, uint8_t* pkind) {
    *pkind = 0;
    if (!strcmp(name, "particles")) return VIAL_LEVEL_PARTICLE;
    if (!strcmp(name, "atoms")) { *pkind = 1; return VIAL_LEVEL_PARTICLE; }
    if (!strcmp(name, "beads")) { *pkind = 2; return VIAL_LEVEL_PARTICLE; }
    if (!strcmp(name, "residues")) return VIAL_LEVEL_RESIDUE;
    if (!strcmp(name, "parts")) return VIAL_LEVEL_PART;
    if (!strcmp(name, "molecules")) return VIAL_LEVEL_MOLECULE;
    if (!strcmp(name, "chains")) return VIAL_LEVEL_CHAIN;
    return VIAL_LEVEL_NONE;
}

static bool k_select(vial_kctx_t* k, const char* name, const bargs_t* ba, vial_val_t* out) {
    if (!k->in->sys) return vial_rt_fail(k, "R0103", "%s needs a system", name);
    uint8_t pkind;
    const uint8_t level = level_of(name, &pkind);
    if (level == VIAL_LEVEL_NONE) {
        // A predefined set
        vial_sets_t* s = vial_select_predefined(k->in->sys, str_from_cstr(name), k->alloc);
        if (!s) return vial_rt_fail(k, "R0103", "%s is not a set", name);
        *out = vial_val_sets(s);
        return true;
    }
    vial_select_t q = {.level = level, .particle_kind = pkind};
    for (uint32_t i = 0; i < ba->n; ++i) {
        const vial_val_t* v = ba->v[i];
        if (!v) continue;
        if (!ba->names[i].len) {
            if (v->dtype != VIAL_DT_SET) return vial_rt_fail(k, "R0103", "%s takes sets as positional arguments here", name);
            q.within = v->sets;
            continue;
        }
        const uint8_t key = vial_select_key(ba->names[i]);
        if (key == VIAL_KEY_SCOPE) {
            if (v->dtype == VIAL_DT_SET) q.scope = v->sets->single ? v->sets : vial_sets_merge(v->sets, k->alloc);
            continue;
        }
        if (key == VIAL_KEY_NONE) return vial_rt_fail(k, "R0103", "%s takes no keyword %.*s", name, (int)ba->names[i].len, ba->names[i].ptr);
        vial_filter_t f = {.key = key};
        const size_t cnt = vial_val_count(v);
        for (size_t e = 0; e < cnt; ++e) {
            vial_match_t m = {0};
            if (v->dtype == VIAL_DT_STR) { m.kind = 0; m.s = ((const str_t*)v->data)[e]; }
            else if (v->dtype == VIAL_DT_RANGE) { m.kind = 2; m.lo = (int64_t)((const double*)v->data)[0]; m.hi = (int64_t)((const double*)v->data)[1]; }
            else if (v->dtype == VIAL_DT_INT) { m.kind = 1; m.lo = ((const int32_t*)v->data)[e]; }
            else return vial_rt_fail(k, "R0103", "a keyword's value here is a string, an int or a range");
            md_array_push(f.any, m, k->alloc);
        }
        md_array_push(q.filters, f, k->alloc);
    }
    vial_sets_t* s = vial_select(k->in->sys, &q, k->alloc);
    if (!s) return vial_rt_fail(k, "R0103", "%s selects nothing it can", name);
    vial_val_t o = vial_val_sets(s);
    // A selector with a positional set that is one set gives one set when the check typed it so
    if (k->type && k->type->ndim == 0 && !s->single) {
        vial_sets_t* m = vial_sets_merge(s, k->alloc);
        m->single = true;
        o = vial_val_sets(m);
    }
    *out = o;
    return true;
}

// ### SERIES OF FRAMES ###

static bool frames_in(vial_kctx_t* k, const vial_val_t* x, const char* name) {
    if (!x || !(x->flags & VIAL_VAL_FRAMES)) return vial_rt_fail(k, "R0103", "%s reads a series of frames", name);
    return true;
}

// Counts of the values in bins over range; values outside it are left out. density=True divides by count and width.
static bool k_histogram(vial_kctx_t* k, const vial_val_t* x, const vial_val_t* bins, const vial_val_t* range, const vial_val_t* density, vial_val_t* out) {
    if (!bins) return vial_rt_fail(k, "R0103", "a missing argument");
    if (!frames_in(k, x, "histogram")) return false;
    const double bd = num_at(bins, 0);
    if (!(bd >= 1.0) || !vial_k_room(k, bd > 1e18 ? SIZE_MAX : (size_t)bd)) return bd >= 1.0 ? false : vial_rt_fail(k, "R0103", "bins= is at least 1");
    const uint32_t nb = (uint32_t)bd;
    if (!range || range->dtype != VIAL_DT_TUPLE || range->num_elems != 2) return vial_rt_fail(k, "R0103", "range= is (lo, hi)");
    const double lo = num_at(&range->elems[0], 0), hi = num_at(&range->elems[1], 0);
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, 1, &nb, &o)) return false;
    double* c = md_alloc(k->alloc, (nb ? nb : 1) * sizeof(double));
    memset(c, 0, nb * sizeof(double));
    const size_t n = vial_val_count(x);
    size_t inside = 0;
    const double w = (hi - lo) / nb;
    for (size_t i = 0; i < n; ++i) {
        if (vial_val_is_missing(x, i)) continue;
        const double v = num_at(x, i);
        if (!(v >= lo && v <= hi)) continue;
        uint32_t b = (uint32_t)((v - lo) / w);
        if (b >= nb) b = nb - 1;    // hi itself falls in the last bin, as NumPy's
        c[b] += 1.0;
        inside += 1;
    }
    const bool dens = density && density->dtype == VIAL_DT_BOOL && *(const uint8_t*)density->data;
    for (uint32_t b = 0; b < nb; ++b) ((float*)o.data)[b] = (float)(dens && inside ? c[b] / ((double)inside * w) : c[b]);
    *out = o;
    return true;
}

// A centred window of n frames that shrinks at the ends, so it stays centred
static bool k_smooth(vial_kctx_t* k, const vial_val_t* x, const vial_val_t* nv, vial_val_t* out) {
    if (!nv) return vial_rt_fail(k, "R0103", "a missing argument");
    if (!frames_in(k, x, "smooth")) return false;
    const int64_t n = (int64_t)num_at(nv, 0), h = n / 2;
    const uint32_t nf = x->dims[0];
    const size_t inner = nf ? vial_val_count(x) / nf : 0;
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, x->ndim, x->dims, &o)) return false;
    o.flags = VIAL_VAL_FRAMES;
    for (uint32_t f = 0; f < nf; ++f) {
        int64_t w = h;
        if (w > f) w = f;
        if (w > (int64_t)nf - 1 - f) w = (int64_t)nf - 1 - f;
        for (size_t e = 0; e < inner; ++e) {
            double s = 0.0;
            uint32_t cnt = 0;
            for (int64_t g = (int64_t)f - w; g <= (int64_t)f + w; ++g) {
                const size_t at = (size_t)g * inner + e;
                if (vial_val_is_missing(x, at)) continue;
                s += num_at(x, at);
                cnt += 1;
            }
            if (cnt) ((float*)o.data)[f * inner + e] = (float)(s / cnt);
            else vial_val_set_missing(&o, f * inner + e, k->alloc);
        }
    }
    *out = o;
    return true;
}

// Central differences over time(); one-sided at the ends
static bool k_derivative(vial_kctx_t* k, const vial_val_t* x, vial_val_t* out) {
    if (!frames_in(k, x, "derivative")) return false;
    const uint32_t nf = x->dims[0];
    const size_t inner = nf ? vial_val_count(x) / nf : 0;
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, x->ndim, x->dims, &o)) return false;
    o.flags = VIAL_VAL_FRAMES;
    const double* t = k->in->times;
    for (uint32_t f = 0; f < nf; ++f) {
        const uint32_t a = f > 0 ? f - 1 : 0, b = f + 1 < nf ? f + 1 : f;
        for (size_t e = 0; e < inner; ++e) {
            const size_t ia = (size_t)a * inner + e, ib = (size_t)b * inner + e;
            if (a == b || vial_val_is_missing(x, ia) || vial_val_is_missing(x, ib) || t[b] == t[a]) { vial_val_set_missing(&o, f * inner + e, k->alloc); continue; }
            ((float*)o.data)[f * inner + e] = (float)((num_at(x, ib) - num_at(x, ia)) / (t[b] - t[a]));
        }
    }
    *out = o;
    return true;
}

// The mean removed, normalized to 1 at lag 0; the elements of an array add up
static bool k_autocorrelation(vial_kctx_t* k, const vial_val_t* x, const vial_val_t* lagv, vial_val_t* out) {
    if (!lagv) return vial_rt_fail(k, "R0103", "a missing argument");
    if (!frames_in(k, x, "autocorrelation")) return false;
    const double ld = num_at(lagv, 0);
    if (!(ld >= 0.0) || !vial_k_room(k, ld > 1e18 ? SIZE_MAX : (size_t)ld + 1)) return ld >= 0.0 ? false : vial_rt_fail(k, "R0103", "max_lag is at least 0");
    const uint32_t nl = (uint32_t)ld + 1;
    const uint32_t nf = x->dims[0];
    const size_t inner = nf ? vial_val_count(x) / nf : 0;
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, 1, &nl, &o)) return false;
    double* mean = md_alloc(k->alloc, (inner ? inner : 1) * sizeof(double));
    for (size_t e = 0; e < inner; ++e) {
        double s = 0.0;
        uint32_t c = 0;
        for (uint32_t f = 0; f < nf; ++f) if (!vial_val_is_missing(x, (size_t)f * inner + e)) { s += num_at(x, (size_t)f * inner + e); c++; }
        mean[e] = c ? s / c : 0.0;
    }
    double c0 = 0.0;
    for (uint32_t lag = 0; lag < nl; ++lag) {
        double s = 0.0;
        uint32_t c = 0;
        for (uint32_t f = 0; f + lag < nf; ++f) {
            for (size_t e = 0; e < inner; ++e) {
                const size_t ia = (size_t)f * inner + e, ib = (size_t)(f + lag) * inner + e;
                if (vial_val_is_missing(x, ia) || vial_val_is_missing(x, ib)) continue;
                s += (num_at(x, ia) - mean[e]) * (num_at(x, ib) - mean[e]);
                c++;
            }
        }
        const double v = c ? s / c : 0.0;
        if (lag == 0) c0 = v;
        if (!c || c0 == 0.0) vial_val_set_missing(&o, lag, k->alloc);
        else ((float*)o.data)[lag] = (float)(v / c0);
    }
    *out = o;
    return true;
}

// ### DISPATCH ###

bool vial_k_builtin(vial_kctx_t* k, const char* name, const vial_val_t* const* args, const str_t* names, uint32_t n, vial_val_t* out) {
    bargs_t ba = {.v = args, .names = names, .n = n};
    // The first letter first: most names differ there, and a mapped call dispatches once per element
    #define IS(s) (name[0] == (s)[0] && strcmp(name, s) == 0)
    if (IS("positions")) return vial_k_positions(k, arg(&ba, "s", 0), out);
    if (IS("com")) return k_com(k, arg(&ba, "s", 0), arg(&ba, "weights", 1), out);
    if (IS("displacement")) return k_displacement(k, arg(&ba, "a", 0), arg(&ba, "b", 1), out);
    if (IS("angle")) return k_angle(k, &ba, 3, out);
    if (IS("dihedral")) return k_angle(k, &ba, 4, out);
    if (IS("closest_pair") || IS("farthest_pair")) {
        const vial_val_t* a[3] = {arg(&ba, "p", 0), arg(&ba, "q", 1), arg(&ba, "distinct", 2)};
        if (!a[0] || !a[1]) return vial_rt_fail(k, "R0103", "%s takes p and q", name);
        k->map_index = -1;
        return mapped(k, a, (const uint8_t[]){VX_KIND_POINTS, VX_KIND_POINTS, VX_KIND_WHOLE}, 3, IS("closest_pair") ? core_pair_c : core_pair_f, out);
    }
    if (IS("fit")) {
        const vial_val_t* a[4] = {arg(&ba, "p", 0), arg(&ba, "ref", 1), arg(&ba, "on", 2), arg(&ba, "weights", 3)};
        return mapped(k, a, (const uint8_t[]){VX_KIND_POINTS, VX_KIND_POINTS, VX_KIND_WHOLE, VX_KIND_WHOLE}, 4, core_fit, out);
    }
    if (IS("within")) return k_within(k, arg(&ba, "r", 0), arg(&ba, "of", 1), out);
    if (IS("eigh")) return k_eigh(k, arg(&ba, "m", 0), out);
    if (IS("trace") || IS("det") || IS("inv") || IS("solve")) return k_linalg(k, name, arg(&ba, "m", 0), arg(&ba, "b", 1), out);
    if (IS("covariance")) {
        const vial_val_t* x = arg(&ba, "x", 0);
        const vial_val_t* axis = arg(&ba, "axis", 1);
        if (x && (x->flags & VIAL_VAL_FRAMES) && axis && axis->dtype == VIAL_DT_STR) {
            // Over frames: the frames are the samples
            if (x->ndim != 2) return vial_rt_fail(k, "R0103", "covariance over frames takes one array of values per frame");
            vial_val_t plain = *x;
            plain.flags &= ~VIAL_VAL_FRAMES;
            return k_covariance(k, &plain, &(vial_val_t){.dtype = VIAL_DT_INT, .data = (int32_t[]){0}}, out);
        }
        return k_covariance(k, x, axis, out);
    }
    if (IS("histogram")) return k_histogram(k, arg(&ba, "x", 0), arg(&ba, "bins", 2), arg(&ba, "range", 3), arg(&ba, "density", 4), out);
    if (IS("smooth")) return k_smooth(k, arg(&ba, "x", 0), arg(&ba, "n", 1), out);
    if (IS("derivative")) return k_derivative(k, arg(&ba, "x", 0), out);
    if (IS("autocorrelation")) return k_autocorrelation(k, arg(&ba, "x", 0), arg(&ba, "max_lag", 1), out);
    if (IS("x") || IS("y") || IS("z")) return k_component(k, args[0], IS("x") ? 0 : IS("y") ? 1 : 2, out);
    if (IS("T")) return k_transpose(k, args[0], out);
    if (IS("reshape")) return k_reshape(k, arg(&ba, "a", 0), out);
    if (IS("flatten")) {
        const vial_val_t* a = arg(&ba, "a", 0);
        vial_val_t o = *a;
        o.ndim = 1;
        o.dims[0] = (uint32_t)vial_val_count(a);
        o.origin = NULL;
        *out = o;
        return true;
    }
    if (IS("linspace")) return k_linspace(k, arg(&ba, "lo", 0), arg(&ba, "hi", 1), arg(&ba, "n", 2), out);
    if (IS("arange")) return k_arange(k, arg(&ba, "lo", 0), arg(&ba, "hi", 1), arg(&ba, "step", 2), out);
    if (IS("eye")) return k_eye(k, arg(&ba, "n", 0), out);
    if (IS("spread")) return k_spread(k, arg(&ba, "values", 0), arg(&ba, "groups", 1), arg(&ba, "combine", 2), out);
    if (IS("mass") || IS("charge") || IS("name") || IS("element")) return k_property(k, name, arg(&ba, "s", 0), out);
    if (IS("norm") || IS("normalize")) return k_norm(k, arg(&ba, "v", 0), IS("normalize"), out);
    if (IS("dot") || IS("cross")) return k_dot_cross(k, arg(&ba, "u", 0), arg(&ba, "v", 1), IS("cross"), out);
    if (IS("sum") || IS("mean") || IS("std") || IS("var") || IS("any") || IS("all") || IS("argmin") || IS("argmax") ||
        ((IS("min") || IS("max")) && !arg(&ba, "b", 1))) {
        const red_t r = IS("sum") ? R_SUM : IS("mean") ? R_MEAN : IS("std") ? R_STD : IS("var") ? R_VAR : IS("any") ? R_ANY :
                        IS("all") ? R_ALL : IS("argmin") ? R_ARGMIN : IS("argmax") ? R_ARGMAX : IS("min") ? R_MIN : R_MAX;
        if (arg(&ba, "by", 99)) return k_reduce_by(k, r, arg(&ba, "a", 0), arg(&ba, "by", 99), out);
        return k_reduce(k, r, arg(&ba, "a", 0), arg(&ba, "axis", 1), out);
    }
    if (IS("min") || IS("max")) {
        const vial_val_t* xs[2] = {arg(&ba, "a", 0), arg(&ba, "b", 1)};
        return k_elementwise(k, name, xs, 2, out);
    }
    if (IS("abs") || IS("floor") || IS("ceil") || IS("round") || IS("sqrt") || IS("exp") || IS("log") || IS("log10") ||
        IS("sin") || IS("cos") || IS("tan") || IS("asin") || IS("acos") || IS("atan")) {
        return k_math1(k, name, arg(&ba, "x", 0), out);
    }
    if (IS("atan2")) { const vial_val_t* xs[2] = {arg(&ba, "y", 0), arg(&ba, "x", 1)}; return k_elementwise(k, name, xs, 2, out); }
    if (IS("clip")) { const vial_val_t* xs[3] = {arg(&ba, "x", 0), arg(&ba, "lo", 1), arg(&ba, "hi", 2)}; return k_elementwise(k, name, xs, 3, out); }
    if (IS("where")) {
        const vial_val_t* xs[3] = {arg(&ba, "c", 0), arg(&ba, "a", 1), arg(&ba, "b", 2)};
        if (!xs[0] || !xs[1] || !xs[2]) return vial_rt_fail(k, "R0103", "where takes three arguments");
        return k_elementwise(k, name, xs, 3, out);
    }
    if (IS("unitless") || IS("float") || IS("int")) return vial_k_convert(k, arg(&ba, "x", 0), k->type->dtype, out);
    if (IS("len")) {
        const vial_val_t* a = arg(&ba, "a", 0);
        const int32_t v = a->dtype == VIAL_DT_SET && a->ndim == 0 ? (int32_t)vial_sets_size(a->sets, 0) : a->ndim ? (int32_t)a->dims[0] : 1;
        *out = vial_val_int(v, k->alloc);
        return true;
    }
    if (IS("count")) {
        const vial_val_t* s = arg(&ba, "s", 0);
        if (!s || s->dtype != VIAL_DT_SET) return vial_rt_fail(k, "R0103", "count takes sets");
        vial_val_t o; if (!vial_k_make(k, VIAL_DT_INT, s->ndim, s->dims, &o)) return false;
        const size_t c = s->ndim ? vial_val_count(s) : 1;
        for (size_t i = 0; i < c; ++i) ((int32_t*)o.data)[i] = (int32_t)vial_sets_size(s->sets, (uint32_t)i);
        *out = o;
        return true;
    }
    if (IS("shown")) { *out = *arg(&ba, "x", 0); return true; }
    if (IS("marker") || IS("line") || IS("arrow") || IS("sphere") || IS("arc") || IS("disk") || IS("label") || IS("highlight")) {
        // Drawings are values for hover; what they draw is made when hover asks. One per element of the broadcast
        // points, as the type says.
        vial_val_t o = {.dtype = VIAL_DT_DRAWING, .ndim = k->type->ndim};
        for (uint8_t d = 0; d < o.ndim; ++d) {
            if (k->type->dims[d].kind == VIAL_EXT_STATIC) { o.dims[d] = k->type->dims[d].n; continue; }
            // An extent known only now: the arguments' leading dimension of that place
            o.dims[d] = 1;
            for (uint32_t i = 0; i < n; ++i) {
                const vial_val_t* v = args[i];
                if (!v || v->ndim < 1) continue;
                const int lead = (v->dtype == VIAL_DT_SET ? v->ndim : v->ndim - 1) - (o.ndim - d);
                if (lead >= 0 && v->dims[lead] > o.dims[d]) o.dims[d] = v->dims[lead];
            }
        }
        *out = o;
        return true;
    }
    if (IS("range")) return k_range(k, &ba, out);
    if (IS("zeros")) return k_fill(k, 0.0, out);
    if (IS("ones")) return k_fill(k, 1.0, out);
    if (IS("full")) { const vial_val_t* v = arg(&ba, "value", 1); return k_fill(k, v ? num_at(v, 0) : 0.0, out); }
    if (IS("zeros_like")) return k_like(k, arg(&ba, "x", 0), 0.0, out);
    if (IS("ones_like")) return k_like(k, arg(&ba, "x", 0), 1.0, out);
    if (IS("concatenate")) return k_concatenate(k, arg(&ba, "arrays", 0), arg(&ba, "axis", 1), out);
    if (IS("missing_where")) return k_missing_where(k, arg(&ba, "c", 0), arg(&ba, "x", 1), out);
    if (IS("valid")) return k_valid(k, arg(&ba, "x", 0), out);
    if (IS("filled")) return k_filled(k, arg(&ba, "x", 0), arg(&ba, "v", 1), out);
    if (IS("diff")) return k_diff(k, arg(&ba, "a", 0), arg(&ba, "axis", 1), out);
    if (IS("frame")) {
        if (!k->in->frame) return vial_rt_fail(k, "R0103", "frame() needs a frame");
        *out = vial_val_int((int32_t)k->in->frame->index, k->alloc);
        return true;
    }
    if (IS("time")) {
        if (!k->in->frame) return vial_rt_fail(k, "R0103", "time() needs a frame");
        *out = vial_val_float64(k->in->frame->time, k->alloc);
        return true;
    }
    if (IS("shape")) {
        const vial_val_t* a = arg(&ba, "a", 0);
        vial_val_t* el = md_alloc(k->alloc, (a->ndim ? a->ndim : 1) * sizeof(vial_val_t));
        for (uint8_t d = 0; d < a->ndim; ++d) el[d] = vial_val_int((int32_t)a->dims[d], k->alloc);
        *out = (vial_val_t){.dtype = VIAL_DT_TUPLE, .elems = el, .num_elems = a->ndim};
        return true;
    }
    if (IS("particles") && n == 1 && !names[0].len && args[0]->dtype != VIAL_DT_SET) {
        // The particles an array runs over
        if (!args[0]->origin) return vial_rt_fail(k, "R0103", "the particles of an array not computed from positions");
        vial_sets_t* s = vial_sets_pick(args[0]->origin, (uint32_t[]){0}, 1, k->alloc);
        s->single = true;
        *out = vial_val_sets(s);
        return true;
    }
    {
        uint8_t pk;
        if (level_of(name, &pk) != VIAL_LEVEL_NONE || IS("protein") || IS("nucleic") || IS("water") || IS("ion") || IS("lipid") ||
            IS("backbone") || IS("sidechain") || IS("system")) {
            return k_select(k, name, &ba, out);
        }
    }
    #undef IS
    return vial_rt_fail(k, "R0103", "%s is not evaluated yet", name);
}

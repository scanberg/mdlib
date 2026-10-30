// The interpreter's built-ins: geometry, properties, reductions, math and the array functions. Each takes its
// arguments by parameter, as the check bound them. Built with floating-point contraction off, like kernels.c.

#include "interp.h"
#include "infer_internal.h"
#include "ast.h"

#include <md_system.h>
#include <md_util.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_vec_math.h>

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
        if (str_eq(a->names[i], str_from_cstr(name))) return a->v[i];
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

static vial_val_t make(vial_kctx_t* k, uint8_t dtype, uint8_t nd, const uint32_t* dims) { return vial_val_make(dtype, nd, dims, k->alloc); }

static bool str_is(const vial_val_t* v, const char* s) {
    return v && v->dtype == VIAL_DT_STR && v->ndim == 0 && str_eq(*(const str_t*)v->data, str_from_cstr(s));
}

// ### POSITIONS ###

static const float* frame_xyz(vial_kctx_t* k) {
    if (!k->in->frame || !k->in->frame->xyz) { vial_rt_fail(k, "R0103", "positions need a frame"); return NULL; }
    return k->in->frame->xyz;
}

bool vial_k_positions(vial_kctx_t* k, const vial_val_t* s, vial_val_t* out) {
    if (s->dtype != VIAL_DT_SET) { *out = *s; return true; }   // coordinates are themselves
    const float* xyz = frame_xyz(k);
    if (!xyz) return false;
    const size_t na = k->in->frame->num_atoms;
    const vial_sets_t* set = s->sets;
    if (s->ndim == 0) {
        const uint32_t n = vial_sets_size(set, 0);
        vial_val_t o = make(k, VIAL_DT_FLOAT, 2, (uint32_t[]){n, 3});
        const uint32_t* idx = vial_sets_items(set, 0);
        for (uint32_t i = 0; i < n; ++i) {
            if (idx[i] >= na) return vial_rt_fail(k, "R0103", "particle %u is beyond the frame's %zu", idx[i], na);
            memcpy((float*)o.data + 3 * i, xyz + 3 * (size_t)idx[i], 3 * sizeof(float));
        }
        o.origin = set;
        *out = o;
        return true;
    }
    // An array of sets of one particle each: one position per set
    const size_t cnt = vial_val_count(s);
    uint32_t dims[VIAL_MAX_DIMS + 1];
    memcpy(dims, s->dims, s->ndim * sizeof(uint32_t));
    dims[s->ndim] = 3;
    vial_val_t o = make(k, VIAL_DT_FLOAT, s->ndim + 1, dims);
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
    vial_val_t o = make(k, VIAL_DT_FLOAT, nd + 1, odims);
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
        vec3_t dx = {{B[0] - A[0], B[1] - A[1], B[2] - A[2]}};
        md_util_min_image_vec3(&dx, 1, &cell);
        float* O = (float*)o.data + 3 * i;
        O[0] = dx.x; O[1] = dx.y; O[2] = dx.z;
        bool miss = false;
        for (int c = 0; c < 3; ++c) miss = miss || vial_val_is_missing(&pa, 3 * ia + c) || vial_val_is_missing(&pb, 3 * ib + c);
        if (miss) for (int c = 0; c < 3; ++c) vial_val_set_missing(&o, 3 * i + c, k->alloc);
    }
    if (pb.origin && nd == lb.ndim) o.origin = pb.origin;
    *out = o;
    return true;
}

// The centre of mass of each set or group of points; periodic-aware
static bool k_com(vial_kctx_t* k, const vial_val_t* s, const vial_val_t* weights, vial_val_t* out) {
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
        vial_val_t o = make(k, VIAL_DT_FLOAT, s->ndim + 1, dims);
        for (size_t i = 0; i < cnt; ++i) {
            const uint32_t n = vial_sets_size(s->sets, (uint32_t)i);
            const int32_t* idx = (const int32_t*)vial_sets_items(s->sets, (uint32_t)i);
            if (!n) { for (int c = 0; c < 3; ++c) vial_val_set_missing(&o, 3 * i + c, k->alloc); continue; }
            float* w = NULL;
            if (by_mass && sys) {
                w = md_alloc(k->alloc, n * sizeof(float));
                for (uint32_t q = 0; q < n; ++q) w[q] = md_atom_mass(&sys->atom, (size_t)idx[q]);
                // The weights are indexed as the positions are: by particle
                float* wf = md_alloc(k->alloc, k->in->frame->num_atoms * sizeof(float));
                for (uint32_t q = 0; q < n; ++q) wf[idx[q]] = w[q];
                w = wf;
            }
            const vec3_t c = md_util_com_compute((const vec3_t*)xyz, w, idx, n, &cell);
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
    vial_val_t o = make(k, VIAL_DT_FLOAT, p.ndim - 1, p.dims);
    o.dims[p.ndim - 2] = 3;
    float* w = NULL;
    if (by_mass) {
        if (!p.origin || !sys || vial_sets_size(p.origin, 0) != n) return vial_rt_fail(k, "R0103", "com by mass of points whose particles are not known; write weights=None");
        w = md_alloc(k->alloc, n * sizeof(float));
        for (uint32_t q = 0; q < n; ++q) w[q] = md_atom_mass(&sys->atom, vial_sets_items(p.origin, 0)[q]);
    }
    for (size_t g = 0; g < groups; ++g) {
        float* O = (float*)o.data + 3 * g;
        if (!n) { for (int c = 0; c < 3; ++c) vial_val_set_missing(&o, 3 * g + c, k->alloc); continue; }
        const vec3_t c = md_util_com_compute((const vec3_t*)((const float*)p.data + 3 * (size_t)n * g), w, NULL, n, &cell);
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
    const double n1[3] = {(double)b1.y * b2.z - (double)b1.z * b2.y, (double)b1.z * b2.x - (double)b1.x * b2.z, (double)b1.x * b2.y - (double)b1.y * b2.x};
    const double n2[3] = {(double)b2.y * b3.z - (double)b2.z * b3.y, (double)b2.z * b3.x - (double)b2.x * b3.z, (double)b2.x * b3.y - (double)b2.y * b3.x};
    const double lb2 = sqrt((double)b2.x * b2.x + (double)b2.y * b2.y + (double)b2.z * b2.z);
    const double m1[3] = {n1[1] * b2.z / lb2 - n1[2] * b2.y / lb2, n1[2] * b2.x / lb2 - n1[0] * b2.z / lb2, n1[0] * b2.y / lb2 - n1[1] * b2.x / lb2};
    const double x = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
    const double y = m1[0] * n2[0] + m1[1] * n2[1] + m1[2] * n2[2];
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
        vial_val_t o = make(k, VIAL_DT_FLOAT, m.ndim - 2, m.dims);
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
    vial_val_t o = make(k, VIAL_DT_FLOAT, nd, dims);
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

// The closest or farthest pair of points of p and q under the minimum image: (distance, i, j)
static bool k_pair(vial_kctx_t* k, const vial_val_t* pv, const vial_val_t* qv, bool closest, vial_val_t* out) {
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
        // Every pair, in order: the first pair at the extreme wins ties. The search lowering replaces this by a
        // spatial search that gives the same bits.
        float best = 0.0f;
        int32_t bi = -1, bj = -1;
        for (uint32_t i = 0; i < n; ++i) {
            const float* a = (const float*)p.data + 3 * i;
            for (uint32_t j = 0; j < m; ++j) {
                const float* b = (const float*)q.data + 3 * j;
                vec3_t dx = {{b[0] - a[0], b[1] - a[1], b[2] - a[2]}};
                md_util_min_image_vec3(&dx, 1, &cell);
                const float d = sqrtf(dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
                if (bi < 0 || (closest ? d < best : d > best)) { best = d; bi = (int32_t)i; bj = (int32_t)j; }
            }
        }
        *(float*)el[0].data = best;
        *(int32_t*)el[1].data = bi;
        *(int32_t*)el[2].data = bj;
    }
    *out = (vial_val_t){.dtype = VIAL_DT_TUPLE, .elems = el, .num_elems = 3};
    return true;
}

// ### PROPERTIES ###

static bool k_property(vial_kctx_t* k, const char* name, const vial_val_t* s, vial_val_t* out) {
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
    vial_val_t o = make(k, str ? VIAL_DT_STR : VIAL_DT_FLOAT, 1, &n);
    const md_attribute_t* charges = !strcmp(name, "charge") ? md_attributes_find(&sys->attributes, STR_LIT("atom/charge")) : NULL;
    const float* cq = charges ? md_attribute_view(charges, MD_ATTRIBUTE_TYPE_F32, 1, 1) : NULL;
    for (uint32_t i = 0; i < n; ++i) {
        const size_t a = idx[i];
        if (!strcmp(name, "mass"))         ((float*)o.data)[i] = md_atom_mass(&sys->atom, a);
        else if (!strcmp(name, "charge"))  ((float*)o.data)[i] = cq ? cq[a] : (float)md_atom_formal_charge(&sys->atom, a);
        else if (!strcmp(name, "name"))    ((str_t*)o.data)[i] = md_atom_name(&sys->atom, a);
        else if (!strcmp(name, "element")) ((str_t*)o.data)[i] = md_util_element_symbol(md_atom_atomic_number(&sys->atom, a));
        else return vial_rt_fail(k, "R0103", "%s is not evaluated yet", name);
    }
    o.origin = s->sets;
    *out = o;
    return true;
}

// ### REDUCTIONS ###

// Which dimensions a reduction takes: axis= as an int, a tuple of ints, or none for all; never the frame here
static bool axes_of(vial_kctx_t* k, const vial_val_t* axis, uint8_t ndim, bool* red) {
    for (uint8_t d = 0; d < ndim; ++d) red[d] = axis == NULL || axis->dtype == VIAL_DT_NONE;
    if (!axis || axis->dtype == VIAL_DT_NONE) return true;
    const vial_val_t* items = axis->dtype == VIAL_DT_TUPLE ? axis->elems : axis;
    const uint32_t n = axis->dtype == VIAL_DT_TUPLE ? axis->num_elems : 1;
    for (uint32_t i = 0; i < n; ++i) {
        if (items[i].dtype == VIAL_DT_STR) return vial_rt_fail(k, "R0103", "axis='frame' is a trajectory reduction");
        if (items[i].dtype != VIAL_DT_INT) return vial_rt_fail(k, "R0103", "axis= takes ints");
        int64_t a = *(const int32_t*)items[i].data;
        if (a < 0) a += ndim;
        if (a < 0 || a >= ndim) return vial_rt_fail(k, "R0103", "axis %d of %u dimensions", *(const int32_t*)items[i].data, ndim);
        red[a] = true;
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
    if (!axes_of(k, axis, a->ndim, rd)) return false;
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
    vial_val_t o = make(k, odt, ond, odims);
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

// ### ELEMENTWISE MATH ###

static bool k_math1(vial_kctx_t* k, const char* name, const vial_val_t* a, vial_val_t* out) {
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
        vial_val_t o = make(k, VIAL_DT_INT, x.ndim, x.dims);
        for (size_t i = 0; i < n; ++i) {
            const int32_t v = ((const int32_t*)x.data)[i];
            if (!strcmp(name, "abs")) { if (v == INT32_MIN) vial_val_set_missing(&o, i, k->alloc); else ((int32_t*)o.data)[i] = v < 0 ? -v : v; }
            else ((int32_t*)o.data)[i] = v;    // floor, ceil and round of ints
        }
        o.missing = o.missing ? o.missing : x.missing;
        *out = o;
        return true;
    }
    const vial_val_t x = vial_val_as(a, odt, k->alloc);
    vial_val_t o = make(k, odt, x.ndim, x.dims);
    for (size_t i = 0; i < n; ++i) {
        const double v = num_at(&x, i) * scale;
        double r;
        if (odt == VIAL_DT_FLOAT) {
            const float f = (float)v;
            float rf;
            if (!strcmp(name, "abs")) rf = fabsf(f);
            else if (!strcmp(name, "floor")) rf = floorf(f);
            else if (!strcmp(name, "ceil")) rf = ceilf(f);
            else if (!strcmp(name, "round")) rf = nearbyintf(f);
            else if (!strcmp(name, "sqrt")) rf = sqrtf(f);
            else if (!strcmp(name, "exp")) rf = expf(f);
            else if (!strcmp(name, "log")) rf = logf(f);
            else if (!strcmp(name, "log10")) rf = log10f(f);
            else if (!strcmp(name, "sin")) rf = sinf(f);
            else if (!strcmp(name, "cos")) rf = cosf(f);
            else if (!strcmp(name, "tan")) rf = tanf(f);
            else if (!strcmp(name, "asin")) rf = asinf(f) * (float)DEG_PER_RAD;
            else if (!strcmp(name, "acos")) rf = acosf(f) * (float)DEG_PER_RAD;
            else if (!strcmp(name, "atan")) rf = atanf(f) * (float)DEG_PER_RAD;
            else return vial_rt_fail(k, "R0103", "%s is not evaluated yet", name);
            ((float*)o.data)[i] = rf;
            continue;
        }
        if (!strcmp(name, "abs")) r = fabs(v);
        else if (!strcmp(name, "floor")) r = floor(v);
        else if (!strcmp(name, "ceil")) r = ceil(v);
        else if (!strcmp(name, "round")) r = nearbyint(v);
        else if (!strcmp(name, "sqrt")) r = sqrt(v);
        else if (!strcmp(name, "exp")) r = exp(v);
        else if (!strcmp(name, "log")) r = log(v);
        else if (!strcmp(name, "log10")) r = log10(v);
        else if (!strcmp(name, "sin")) r = sin(v);
        else if (!strcmp(name, "cos")) r = cos(v);
        else if (!strcmp(name, "tan")) r = tan(v);
        else if (!strcmp(name, "asin")) r = asin(v) * DEG_PER_RAD;
        else if (!strcmp(name, "acos")) r = acos(v) * DEG_PER_RAD;
        else if (!strcmp(name, "atan")) r = atan(v) * DEG_PER_RAD;
        else return vial_rt_fail(k, "R0103", "%s is not evaluated yet", name);
        ((double*)o.data)[i] = r;
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
    vial_val_t o = make(k, odt, nd, dims);
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
    const vial_val_t x = vial_val_as(a, VIAL_DT_FLOAT, k->alloc);
    if (!x.ndim) return vial_rt_fail(k, "R0103", "norm of a scalar");
    const uint32_t m = x.dims[x.ndim - 1];
    const size_t g = m ? vial_val_count(&x) / m : 0;
    vial_val_t o = normalize ? make(k, VIAL_DT_FLOAT, x.ndim, x.dims) : make(k, VIAL_DT_FLOAT, x.ndim - 1, x.dims);
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
    vial_val_t o = make(k, VIAL_DT_FLOAT, cross ? nd + 1 : nd, dims);
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
    vial_val_t o = make(k, t->dtype, t->ndim, dims);
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
    vial_val_t o = make(k, k->type->dtype, a->ndim, a->dims);
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
    const uint32_t n = (uint32_t)c;
    vial_val_t o = make(k, VIAL_DT_INT, 1, &n);
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
    vial_val_t o = make(k, dt, nd, dims);
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
    vial_val_t o = make(k, VIAL_DT_BOOL, x->ndim, x->dims);
    const size_t n = vial_val_count(x);
    for (size_t i = 0; i < n; ++i) ((uint8_t*)o.data)[i] = !vial_val_is_missing(x, i);
    *out = o;
    return true;
}

static bool k_filled(vial_kctx_t* k, const vial_val_t* x, const vial_val_t* v, vial_val_t* out) {
    const uint8_t dt = k->type->dtype;
    vial_val_t o = vial_val_as(x, dt, k->alloc);
    vial_val_t c = make(k, dt, o.ndim, o.dims);
    const size_t n = vial_val_count(&o), es = vial_dtype_size(dt);
    memcpy(c.data, o.data, n * es);
    const vial_val_t fv = vial_val_as(v, dt, k->alloc);
    for (size_t i = 0; i < n; ++i) if (vial_val_is_missing(x, i)) memcpy((uint8_t*)c.data + i * es, fv.data, es);
    *out = c;
    return true;
}

static bool k_diff(vial_kctx_t* k, const vial_val_t* a, const vial_val_t* axis, vial_val_t* out) {
    const uint8_t dt = k->type->dtype;
    const vial_val_t x = vial_val_as(a, dt, k->alloc);
    int64_t ax = axis && axis->dtype == VIAL_DT_INT ? *(const int32_t*)axis->data : 0;
    if (ax < 0) ax += x.ndim;
    if (ax < 0 || ax >= x.ndim) return vial_rt_fail(k, "R0103", "axis out of range");
    uint32_t dims[VIAL_MAX_DIMS];
    memcpy(dims, x.dims, x.ndim * sizeof(uint32_t));
    if (!dims[ax]) return vial_rt_fail(k, "R0103", "diff of an empty dimension");
    dims[ax] -= 1;
    vial_val_t o = make(k, dt, x.ndim, dims);
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

// ### DISPATCH ###

bool vial_k_builtin(vial_kctx_t* k, const char* name, const vial_val_t* const* args, const str_t* names, uint32_t n, vial_val_t* out) {
    bargs_t ba = {.v = args, .names = names, .n = n};
    #define IS(s) (strcmp(name, s) == 0)
    if (IS("positions")) return vial_k_positions(k, arg(&ba, "s", 0), out);
    if (IS("com")) return k_com(k, arg(&ba, "s", 0), arg(&ba, "weights", 1), out);
    if (IS("displacement")) return k_displacement(k, arg(&ba, "a", 0), arg(&ba, "b", 1), out);
    if (IS("angle")) return k_angle(k, &ba, 3, out);
    if (IS("dihedral")) return k_angle(k, &ba, 4, out);
    if (IS("closest_pair") || IS("farthest_pair")) return k_pair(k, arg(&ba, "p", 0), arg(&ba, "q", 1), IS("closest_pair"), out);
    if (IS("mass") || IS("charge") || IS("name") || IS("element")) return k_property(k, name, arg(&ba, "s", 0), out);
    if (IS("norm") || IS("normalize")) return k_norm(k, arg(&ba, "v", 0), IS("normalize"), out);
    if (IS("dot") || IS("cross")) return k_dot_cross(k, arg(&ba, "u", 0), arg(&ba, "v", 1), IS("cross"), out);
    if (IS("sum") || IS("mean") || IS("std") || IS("var") || IS("any") || IS("all") || IS("argmin") || IS("argmax") ||
        ((IS("min") || IS("max")) && !arg(&ba, "b", 1))) {
        const red_t r = IS("sum") ? R_SUM : IS("mean") ? R_MEAN : IS("std") ? R_STD : IS("var") ? R_VAR : IS("any") ? R_ANY :
                        IS("all") ? R_ALL : IS("argmin") ? R_ARGMIN : IS("argmax") ? R_ARGMAX : IS("min") ? R_MIN : R_MAX;
        if (arg(&ba, "by", 99)) return vial_rt_fail(k, "R0103", "by= is not evaluated yet");
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
        vial_val_t o = make(k, VIAL_DT_INT, s->ndim, s->dims);
        const size_t c = s->ndim ? vial_val_count(s) : 1;
        for (size_t i = 0; i < c; ++i) ((int32_t*)o.data)[i] = (int32_t)vial_sets_size(s->sets, (uint32_t)i);
        *out = o;
        return true;
    }
    if (IS("shown")) { *out = *arg(&ba, "x", 0); return true; }
    if (IS("marker") || IS("line") || IS("arrow") || IS("sphere") || IS("arc") || IS("disk") || IS("label") || IS("highlight")) {
        // Drawings are values for hover; what they draw is made when hover asks
        *out = (vial_val_t){.dtype = VIAL_DT_DRAWING};
        return true;
    }
    if (IS("range")) return k_range(k, &ba, out);
    if (IS("zeros")) return k_fill(k, 0.0, out);
    if (IS("ones")) return k_fill(k, 1.0, out);
    if (IS("full")) { const vial_val_t* v = arg(&ba, "value", 1); return k_fill(k, v ? num_at(v, 0) : 0.0, out); }
    if (IS("zeros_like")) return k_like(k, arg(&ba, "a", 0), 0.0, out);
    if (IS("ones_like")) return k_like(k, arg(&ba, "a", 0), 1.0, out);
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
    #undef IS
    return vial_rt_fail(k, "R0103", "%s is not evaluated yet", name);
}

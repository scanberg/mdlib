// Volumes: sdf() and density() (volume.h)

#include "volume.h"
#include "builtins.h"

#include <md_system.h>
#include <md_geom.h>
#include <md_mcs.h>
#include <md_unitcell.h>
#include <core/md_spatial_acc.h>
#include <core/md_os.h>
#include <core/md_coord_stream.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_vec_math.h>

#include <float.h>
#include <math.h>
#include <string.h>

bool vial_k_positions(vial_kctx_t* k, const vial_val_t* s, vial_val_t* out);

// ### CORRESPONDENCE ###

// Whether two structures are alike: as many particles, of the same elements, or names where an element is not known
static bool alike_structures(const md_system_t* sys, const uint32_t* a, const uint32_t* b, uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) {
        if (a[i] >= sys->atom.count || b[i] >= sys->atom.count) return false;
        const md_atomic_number_t za = md_atom_atomic_number(&sys->atom, a[i]), zb = md_atom_atomic_number(&sys->atom, b[i]);
        if ((za && zb) ? za != zb : !md_str_eq(md_atom_name(&sys->atom, a[i]), md_atom_name(&sys->atom, b[i]))) return false;
    }
    return true;
}

// The position of atom in an ascending set, or -1
static int pos_in(const uint32_t* items, uint32_t n, int32_t atom) {
    uint32_t lo = 0, hi = n;
    while (lo < hi) { const uint32_t mid = (lo + hi) / 2; if ((int32_t)items[mid] < atom) lo = mid + 1; else hi = mid; }
    return lo < n && (int32_t)items[lo] == atom ? (int)lo : -1;
}

// One structure's mapping onto the first, as positions: col_of[p] for each particle position p of the first is the
// position in the structure that plays its part, or -1. Shared by alike structures.
typedef struct mcs_map_t {
    uint32_t    first_of;   // the structure it was found on
    int32_t*    col_of;     // [n0]
    size_t      pairs, heavy;
    bool        incomplete;
} mcs_map_t;

bool vial_mcs(vial_mcs_t* out, const md_system_t* sys, const vial_sets_t* st, bool by_name, const float* xyz, md_allocator_i* alloc) {
    memset(out, 0, sizeof(*out));
    if (!sys || !st || !st->count) return false;
    const uint32_t k = st->count, n0 = vial_sets_size(st, 0);
    const uint32_t* first = vial_sets_items(st, 0);
    mcs_map_t* maps = md_alloc(alloc, k * sizeof(mcs_map_t));
    uint32_t num_maps = 0;
    uint32_t* map_of = md_alloc(alloc, k * sizeof(uint32_t));
    uint8_t* in_all = md_alloc(alloc, (n0 ? n0 : 1));
    memset(in_all, 1, n0);
    for (uint32_t i = 0; i < k; ++i) {
        const uint32_t n = vial_sets_size(st, i);
        const uint32_t* items = vial_sets_items(st, i);
        // An alike structure already mapped: the same positions
        uint32_t m = UINT32_MAX;
        for (uint32_t q = 0; q < num_maps && m == UINT32_MAX; ++q) {
            const uint32_t f = maps[q].first_of;
            if (vial_sets_size(st, f) == n && alike_structures(sys, vial_sets_items(st, f), items, n)) m = q;
        }
        if (m == UINT32_MAX) {
            const md_mcs_desc_t d = {
                .a = {.sys = sys, .atoms = (const int32_t*)first, .count = n0, .xyz = (const md_float3_t*)xyz},
                .b = {.sys = sys, .atoms = (const int32_t*)items, .count = n, .xyz = (const md_float3_t*)xyz},
                .flags = by_name ? MD_MCS_FLAG_BY_NAME : 0,
            };
            md_mcs_result_t r;
            if (!md_mcs(&r, &d, md_get_heap_allocator())) return false;
            mcs_map_t* mp = &maps[num_maps];
            mp->first_of = i;
            mp->col_of = md_alloc(alloc, (n0 ? n0 : 1) * sizeof(int32_t));
            for (uint32_t p = 0; p < n0; ++p) mp->col_of[p] = -1;
            for (size_t j = 0; j < r.count; ++j) {
                const int pa = pos_in(first, n0, r.a[j]), pb = pos_in(items, n, r.b[j]);
                if (pa >= 0 && pb >= 0) mp->col_of[pa] = pb;
            }
            mp->pairs = r.count;
            mp->heavy = r.a_heavy < r.b_heavy ? r.a_heavy : r.b_heavy;
            mp->incomplete = r.incomplete;
            if (i == 0) out->first_heavy = r.a_heavy;
            md_mcs_result_free(&r);
            m = num_maps++;
            if (i && (!out->least_at || mp->pairs < out->least_common)) { out->least_at = i + 1; out->least_common = mp->pairs; out->least_heavy = mp->heavy; }
            out->incomplete = out->incomplete || mp->incomplete;
        }
        map_of[i] = m;
        for (uint32_t p = 0; p < n0; ++p) if (maps[m].col_of[p] < 0) in_all[p] = 0;
    }
    // The columns: the first's particles every structure has a part for, ascending; then each row
    uint32_t w = 0;
    for (uint32_t p = 0; p < n0; ++p) w += in_all[p];
    int32_t* flat = md_alloc(alloc, ((size_t)k * w + 1) * sizeof(int32_t));
    for (uint32_t i = 0; i < k; ++i) {
        const uint32_t* items = vial_sets_items(st, i);
        const mcs_map_t* mp = &maps[map_of[i]];
        uint32_t c = 0;
        for (uint32_t p = 0; p < n0; ++p) if (in_all[p]) flat[(size_t)i * w + c++] = (int32_t)items[mp->col_of[p]];
    }
    out->rows = vial_sets_particles(flat, k * w, alloc);
    out->k = k;
    out->width = w;
    return true;
}

static const vial_val_t* varg_k(const vial_val_t* const* args, const md_str_t* names, uint32_t n, const char* name, uint32_t pos);

bool vial_k_mcs(vial_kctx_t* k, const vial_val_t* const* args, const md_str_t* names, uint32_t n, const float* xyz, vial_val_t* out) {
    // Every positional argument is structures; match= and least= are named
    const vial_sets_t* parts[64];
    uint32_t np = 0;
    const vial_val_t* by = NULL;
    for (uint32_t i = 0; i < n; ++i) {
        if (names && names[i].len) { if (md_str_eq(names[i], MD_STR_LIT("match"))) by = args[i]; continue; }
        if (!args[i] || args[i]->dtype != VIAL_DT_SET || !args[i]->sets || !args[i]->sets->count) return vial_rt_fail(k, "R0103", "mcs takes sets, or arrays of sets");
        if (np < 64) parts[np++] = args[i]->sets;
    }
    if (!np) return vial_rt_fail(k, "R0103", "mcs takes at least one structure");
    const vial_sets_t* st = np == 1 ? parts[0] : vial_sets_concat(parts, np, k->alloc);
    const bool by_name = by && by->dtype == VIAL_DT_STR && by->data && md_str_eq(*(const md_str_t*)by->data, MD_STR_LIT("name"));
    vial_mcs_t m;
    if (!vial_mcs(&m, k->in->msys, st, by_name, xyz, k->alloc)) return vial_rt_fail(k, "R0103", "mcs: a structure is too large to search");
    *out = vial_val_sets(m.rows);
    out->ndim = 2;
    out->dims[0] = m.k;
    out->dims[1] = m.width;
    return true;
}

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
    const md_str_t*      names;
    uint32_t          n;
} vargs_t;

// By parameter name when the check bound them, else the pos-th positional one; NULL when not given
static const vial_val_t* varg(const vargs_t* a, const char* name, uint32_t pos) {
    bool bound = false;
    for (uint32_t i = 0; i < a->n; ++i) {
        if (!a->names || !a->names[i].len) continue;
        if (md_str_eq(a->names[i], md_str_from_cstr(name))) return a->v[i];
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

static const vial_val_t* varg_k(const vial_val_t* const* args, const md_str_t* names, uint32_t n, const char* name, uint32_t pos) {
    const vargs_t a = {args, names, n};
    return varg(&a, name, pos);
}

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

static int cmp_u32(const void* a, const void* b) {
    const uint32_t x = *(const uint32_t*)a, y = *(const uint32_t*)b;
    return (x > y) - (x < y);
}

// The positions and masses of particles idx[n] of xyz, as v1 weighs a structure
static void gather(md_float4_t* out, const float* xyz, size_t num_atoms, const uint32_t* idx, uint32_t n, const md_system_t* sys) {
    for (uint32_t i = 0; i < n; ++i) {
        const uint32_t p = idx[i];
        const float* q = p < num_atoms ? xyz + 3 * (size_t)p : NULL;
        const float m = sys && p < sys->atom.count ? md_atom_mass(&sys->atom, p) : 1.0f;
        out[i] = q ? md_float4_set(q[0], q[1], q[2], m > 0.0f ? m : 1.0f) : md_float4_set(0, 0, 0, 0);
    }
}

// One kind's structures from its argument: one set, an array of sets alike, or set[k][w] paired by column
static bool kind_from(vial_kctx_t* k, vial_sdf_kind_t* kind, const vial_val_t* v, md_allocator_i* alloc) {
    memset(kind, 0, sizeof(*kind));
    kind->R = md_float3x3_ident();
    if (!v || v->dtype != VIAL_DT_SET || !v->sets || !v->sets->count) return vial_rt_fail(k, "R0103", "sdf takes a set, an array of sets or set[k][w] of structures");
    const vial_sets_t* s = v->sets;
    if (v->ndim == 2) {
        kind->count = v->dims[0];
        kind->n = v->dims[1];
        if ((size_t)kind->count * kind->n != s->count) return vial_rt_fail(k, "R0103", "set[k][w] that does not hold k * w sets");
        kind->atoms = md_alloc(alloc, ((size_t)s->count + 1) * sizeof(uint32_t));
        for (uint32_t i = 0; i < s->count; ++i) {
            if (vial_sets_size(s, i) != 1) return vial_rt_fail(k, "R0103", "the structures of sdf as set[k][w] are of one particle each");
            kind->atoms[i] = vial_sets_items(s, i)[0];
        }
    } else {
        kind->count = s->count;
        kind->n = vial_sets_size(s, 0);
        kind->atoms = md_alloc(alloc, ((size_t)kind->count * kind->n + 1) * sizeof(uint32_t));
        for (uint32_t i = 0; i < s->count; ++i) {
            if (vial_sets_size(s, i) != kind->n) return vial_rt_fail(k, "R0103", "the structures of sdf have different sizes");
            memcpy(kind->atoms + (size_t)i * kind->n, vial_sets_items(s, i), kind->n * sizeof(uint32_t));
        }
    }
    if (!kind->n) return vial_rt_fail(k, "R0103", "the structures of sdf are empty");
    kind->own = md_alloc(alloc, ((size_t)kind->count * kind->n + 1) * sizeof(uint32_t));
    memcpy(kind->own, kind->atoms, (size_t)kind->count * kind->n * sizeof(uint32_t));
    for (uint32_t i = 0; i < kind->count; ++i) qsort(kind->own + (size_t)i * kind->n, kind->n, sizeof(uint32_t), cmp_u32);
    return true;
}

// A kind's reference: ref= (its positions in the kind's order) or the kind's first structure at frame 0, by mass, in
// mutually consistent images about its own centre, as v1 places it
static bool kind_ref(vial_kctx_t* k, vial_sdf_kind_t* kind, const vial_val_t* ref, const vial_frame_t* first, const md_system_t* sys, md_allocator_i* alloc) {
    kind->ref_xyzw = md_alloc(alloc, kind->n * sizeof(md_float4_t));
    if (!is_none(ref)) {
        vial_val_t r = *ref;
        if (r.dtype == VIAL_DT_FLOAT64) r = vial_val_as(ref, VIAL_DT_FLOAT, k->alloc);
        if (r.dtype != VIAL_DT_FLOAT || vial_val_count(&r) != 3 * (size_t)kind->n || !r.data) {
            return vial_rt_fail(k, "R0103", "ref= has the positions of one structure: float[%u][3]", kind->n);
        }
        for (uint32_t i = 0; i < kind->n; ++i) {
            const float* q = (const float*)r.data + 3 * (size_t)i;
            const float m = sys && kind->atoms[i] < sys->atom.count ? md_atom_mass(&sys->atom, kind->atoms[i]) : 1.0f;
            kind->ref_xyzw[i] = md_float4_set(q[0], q[1], q[2], m > 0.0f ? m : 1.0f);
        }
    } else {
        if (!first || !first->xyz) return vial_rt_fail(k, "R0104", "sdf fits onto frame 0, which the run does not have");
        gather(kind->ref_xyzw, first->xyz, first->num_atoms, kind->atoms, kind->n, sys);
    }
    const md_unitcell_t cell = first ? first->cell : md_unitcell_none();
    md_geom_deperiodize_self(kind->ref_xyzw, kind->n, &cell, &kind->ref_com);
    return true;
}

// The structures and their reference: one kind, fitted particle by particle, in the reference's principal axes
static bool sdf_structures(vial_kctx_t* k, vial_volume_t* vol, const vial_val_t* st, const vial_val_t* ref, const vial_frame_t* first,
                           md_allocator_i* alloc) {
    const md_system_t* sys = k->in->msys;
    vol->sys = sys;
    vol->num_kinds = 1;
    vol->kinds = md_alloc(alloc, sizeof(vial_sdf_kind_t));
    vial_sdf_kind_t* kind = &vol->kinds[0];
    if (!kind_from(k, kind, st, alloc)) return false;
    if (!kind_ref(k, kind, ref, first, sys, alloc)) return false;
    vol->num_structures = kind->count;
    vol->max_n = kind->n;
    const md_float3x3_eigen_t eig = md_float3x3_eigen(md_float3x3_covariance_matrix_float4(kind->ref_xyzw, 0, kind->n, kind->ref_com));
    vol->axes = md_float3x3_transpose(eig.vectors);
    return true;
}

bool vial_volume_begin(vial_kctx_t* k, vial_volume_t* vol, const char* name, const vial_val_t* const* args, const md_str_t* names,
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
        if (!three(ext, vol->extent) || !(vol->extent[0] > 0.0 && vol->extent[1] > 0.0 && vol->extent[2] > 0.0)) {
            return vial_rt_fail(k, "R0103", "extent= is a length more than 0, or three");
        }
        if (!sdf_structures(k, vol, st, ref, first, alloc)) return false;
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
    // Its kinds, their structures and references, are read only: shared with vol, which outlives it
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

// Structure s, over the kinds: its kind and its place in it
static const vial_sdf_kind_t* kind_at(const vial_volume_t* vol, uint32_t s, uint32_t* j) {
    for (uint32_t i = 0; i < vol->num_kinds; ++i) {
        if (s < vol->kinds[i].count) { *j = s; return &vol->kinds[i]; }
        s -= vol->kinds[i].count;
    }
    return NULL;
}

const uint32_t* vial_volume_structure(const vial_volume_t* vol, uint32_t s, uint32_t* n) {
    uint32_t j = 0;
    const vial_sdf_kind_t* kind = vol && vol->kind == VIAL_VOLUME_SDF ? kind_at(vol, s, &j) : NULL;
    if (!kind) { *n = 0; return NULL; }
    *n = kind->n;
    return kind->own + (size_t)j * kind->n;
}

// Structure j of a kind fitted at positions xyz: the rotation onto the kind's reference and its centre, in the image
// it arrived in. False when it is too spread over the cell to place.
static bool fit_structure(const vial_sdf_kind_t* kind, uint32_t j, const float* xyz, size_t num_atoms, const md_unitcell_t* cell,
                          const md_system_t* sys, md_float4_t* scratch, md_float3x3_t* R, md_float3_t* com) {
    gather(scratch, xyz, num_atoms, kind->atoms + (size_t)j * kind->n, kind->n, sys);
    const float residual = md_geom_optimal_rotation_pbc(R, com, NULL, kind->ref_xyzw, kind->ref_com, scratch, kind->n, cell, 8, 1.0e-6f);
    // The residual is the margin to an ambiguous choice of images: near half a cell no choice per point represents
    // this copy, so it is left out rather than smeared into the volume in a wrong orientation
    float min_ext = FLT_MAX;
    if (has_cell(cell)) {
        float ext[3];
        md_unitcell_diag_extract_float(ext, cell);
        for (int q = 0; q < 3; ++q) if (ext[q] > 0.0f) min_ext = MD_MIN(min_ext, ext[q]);
    }
    return residual <= 0.25f * min_ext;
}

// The run's frame to the volume's for a structure of a kind fitted by R about com: axes (Rk (R (p - com)) + tk), as
// Q (p - com) + o
static void structure_frame(const vial_volume_t* vol, const vial_sdf_kind_t* kind, md_float3x3_t R, md_float3x3_t* Q, md_float3_t* o) {
    *Q = md_float3x3_mul(vol->axes, md_float3x3_mul(kind->R, R));
    *o = md_float3x3_mul_float3(vol->axes, kind->t);
}

bool vial_volume_pose(const vial_volume_t* vol, uint32_t s, const float* xyz, size_t num_atoms, const md_unitcell_t* cell, md_float4x4_t* out) {
    uint32_t j = 0;
    const vial_sdf_kind_t* kind = vol && vol->kind == VIAL_VOLUME_SDF ? kind_at(vol, s, &j) : NULL;
    if (!kind) return false;
    md_allocator_i* heap = md_get_heap_allocator();
    md_float4_t* scratch = md_alloc(heap, kind->n * sizeof(md_float4_t));
    md_float3x3_t R;
    md_float3_t com;
    const bool ok = fit_structure(kind, j, xyz, num_atoms, cell, vol->sys, scratch, &R, &com);
    md_free(heap, scratch, kind->n * sizeof(md_float4_t));
    if (!ok) return false;
    md_float3x3_t Q;
    md_float3_t o;
    structure_frame(vol, kind, R, &Q, &o);
    const md_float3_t qc = md_float3x3_mul_float3(Q, com);
    *out = md_float4x4_from_float3x3(Q);
    out->elem[3][0] = o.x - qc.x;
    out->elem[3][1] = o.y - qc.y;
    out->elem[3][2] = o.z - qc.z;
    return true;
}

typedef struct sdf_cb_t {
    const vial_volume_t*    vol;
    md_float3x3_t                  Q;          // the run's frame to the volume's
    md_float3_t                  centre;     // where the points come in around
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
        const md_float3_t d = {{x[i] - c->centre.x, y[i] - c->centre.y, z[i] - c->centre.z}};
        const md_float3_t l = md_float3x3_mul_float3(c->Q, d);
        // Grid positions, voxel centres at whole numbers
        const double gx = ((double)l.x - vol->origin[0]) * c->inv_step[0] - 0.5;
        const double gy = ((double)l.y - vol->origin[1]) * c->inv_step[1] - 0.5;
        const double gz = ((double)l.z - vol->origin[2]) * c->inv_step[2] - 0.5;
        if (gx < -1.0 || gy < -1.0 || gz < -1.0 || gx >= vol->dims[0] || gy >= vol->dims[1] || gz >= vol->dims[2]) continue;
        deposit(vol->sum, vol->dims, gx, gy, gz, w, false);
    }
}

bool vial_volume_add(vial_kctx_t* k, vial_volume_t* vol, const vial_val_t* const* args, const md_str_t* names, uint32_t n) {
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
    if (!is_none(wv) && wv->dtype == VIAL_DT_STR) {
        // weights='mass': the particles' masses, so the points must be particles
        const md_system_t* sys = k->in->msys;
        if (!ids || !sys) return vial_rt_fail(k, "R0103", "weights='mass' needs points whose particles are known");
        for (size_t i = 0; i < np; ++i) {
            const double x = ids[i] < sys->atom.count ? md_atom_mass(&sys->atom, ids[i]) : 0.0;
            w[i] = (float)(x * scale);
            total += x;
        }
    } else if (!is_none(wv)) {
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
        md_float4_t* scratch = md_alloc(k->alloc, vol->max_n * sizeof(md_float4_t));
        double reach = 0.0;
        for (int q = 0; q < 3; ++q) reach = MD_MAX(reach, vol->extent[q] + sqrt(vol->basis[q][0] * vol->basis[q][0] + vol->basis[q][1] * vol->basis[q][1] + vol->basis[q][2] * vol->basis[q][2]));
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
        for (uint32_t s = 0; s < vol->num_structures; ++s) {
            uint32_t j = 0;
            const vial_sdf_kind_t* kind = kind_at(vol, s, &j);
            md_float3x3_t R;
            md_float3_t com;
            if (!fit_structure(kind, j, fr->xyz, fr->num_atoms, &cell, k->in->msys, scratch, &R, &com)) continue;
            added += 1;
            if (!np) continue;
            // The volume's box is about cen in the run's frame, Q (p - cen) in the volume's. The points come in around
            // the query's centre, folded into the cell: a lattice vector off cen when the structure sits outside it.
            md_float3x3_t Q;
            md_float3_t o;
            structure_frame(vol, kind, R, &Q, &o);
            const md_float3_t back = md_float3x3_mul_float3(md_float3x3_transpose(Q), o);
            const double cen[3] = {com.x - back.x, com.y - back.y, com.z - back.z};
            double qc[3];
            md_spatial_acc_aabb_query_center(qc, &acc, cen);
            sdf_cb_t c = {.vol = vol, .Q = Q, .centre = {{(float)qc[0], (float)qc[1], (float)qc[2]}},
                          .own = kind->own + (size_t)j * kind->n, .num_own = kind->n, .ids = ids, .w = w};
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
            md_float3_t d = {{(float)(q[0] - c[0]), (float)(q[1] - c[1]), (float)(q[2] - c[2])}};
            if (periodic) md_geom_min_image(d.elem, 0, 1, &cell);
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

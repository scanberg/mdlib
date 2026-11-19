// rdf(): the radial distribution function, added up over the frames (accum.h)

#include "accum.h"

#include <core/md_spatial_acc.h>
#include <core/md_unitcell.h>
#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_unit.h>

#include <math.h>

bool vial_k_positions(vial_kctx_t* k, const vial_val_t* s, vial_val_t* out);

// ### UNITS AND THE RANGE ###

// What one unit of a value of type t is in Å: 1 for a value without a unit, or one not known, which is read as Å
static double length_scale(const vial_type_t* t) {
    if (!t || t->unit_kind != VIAL_UNIT_KNOWN || md_unit_is_none(t->unit)) return 1.0;
    md_unit_t angstrom;
    vial_unit_name(MD_STR_LIT("angstrom"), &angstrom);
    double f;
    return md_unit_conversion_factor(&f, t->unit, angstrom) ? f : 1.0;
}

// The type the check gave argument i of the op, NULL when there is none
static const vial_type_t* arg_type(const vial_kctx_t* k, int i) {
    if (i < 0 || !k->op || !k->in || !k->in->ir || (uint32_t)i >= k->op->num_args) return NULL;
    const vial_ir_t* ir = k->in->ir;
    const vial_ssa_t s = ir->args[k->op->arg_beg + (uint32_t)i];
    return s ? &ir->values[s].type : NULL;
}

bool vial_rdf_range(const vial_val_t* cutoff, double scale, double scale2, double* lo, double* hi) {
    if (!cutoff) return false;
    double l = 0.0, h = 0.0;
    if (cutoff->dtype == VIAL_DT_RANGE && cutoff->data) {
        l = ((const double*)cutoff->data)[0] * scale;
        h = ((const double*)cutoff->data)[1] * scale;
    } else if (cutoff->dtype == VIAL_DT_TUPLE && cutoff->num_elems == 2) {
        const vial_val_t* e = cutoff->elems;
        if (!vial_arg_is_num(&e[0]) || !vial_arg_is_num(&e[1]) || vial_val_count(&e[0]) != 1 || vial_val_count(&e[1]) != 1) return false;
        l = vial_arg_num(&e[0], 0) * scale;
        h = vial_arg_num(&e[1], 0) * scale2;
    } else if (vial_arg_is_num(cutoff) && vial_val_count(cutoff) == 1) {
        h = vial_arg_num(cutoff, 0) * scale;
    } else {
        return false;
    }
    if (!(l >= 0.0 && h > l && isfinite(h))) return false;
    *lo = l;
    *hi = h;
    return true;
}

// ### THE CELL ###

static bool has_cell(const md_unitcell_t* cell) {
    return (md_unitcell_flags(cell) & (MD_UNITCELL_ORTHO | MD_UNITCELL_TRICLINIC)) != 0;
}

static void cross(double out[3], const double a[3], const double b[3]) {
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

// The cell's volume, and its narrowest width across a periodic axis: the distance between the two faces the other two
// basis vectors span. A pair further apart than half of it may be nearer through another image than the minimum image
// finds, so an rdf reaches at most that far.
static double cell_volume_width(const md_unitcell_t* cell, double* width) {
    double A[3][3];
    md_unitcell_A_extract_double(A, cell);
    double c[3];
    cross(c, A[1], A[2]);
    const double vol = fabs(A[0][0] * c[0] + A[0][1] * c[1] + A[0][2] * c[2]);
    double w = INFINITY;
    const md_unitcell_flags_t pbc[3] = {MD_UNITCELL_PBC_X, MD_UNITCELL_PBC_Y, MD_UNITCELL_PBC_Z};
    for (int q = 0; q < 3; ++q) {
        if (!(cell->flags & pbc[q])) continue;
        cross(c, A[(q + 1) % 3], A[(q + 2) % 3]);
        const double area = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
        if (area > 0.0 && vol / area < w) w = vol / area;
    }
    *width = w;
    return vol;
}

// ### THE POINTS ###

typedef struct rdf_points_t {
    const float*    xyz;    // [n][3], Å
    const uint32_t* ids;    // [n]: the particle each point is, NULL when they are not particles
    size_t          n;
} rdf_points_t;

// The points of argument name at this frame, in Å: a set's particles, or points; those missing are left out
static bool points_of(vial_kctx_t* k, const vial_args_t* a, const char* name, uint32_t pos, rdf_points_t* out) {
    const int i = vial_arg_index(a, name, pos);
    const vial_val_t* v = i < 0 ? NULL : a->v[i];
    if (!v) return vial_rt_fail(k, "R0103", "a missing argument");
    vial_val_t p;
    if (!vial_k_positions(k, v, &p)) return false;
    if (p.dtype == VIAL_DT_FLOAT64) p = vial_val_as(&p, VIAL_DT_FLOAT, k->alloc);
    if (p.dtype != VIAL_DT_FLOAT || p.ndim != 2 || p.dims[1] != 3) return vial_rt_fail(k, "R0103", "%s= of `rdf` is one set or points, float[n][3]", name);
    const size_t n = p.dims[0];
    // Which particle each point is: a set's, or what the points are positions of
    const uint32_t* ids = NULL;
    if (v->dtype == VIAL_DT_SET && v->sets && v->sets->count == 1 && vial_sets_size(v->sets, 0) == n) ids = vial_sets_items(v->sets, 0);
    else if (p.origin && p.origin->count == 1 && vial_sets_size(p.origin, 0) == n) ids = vial_sets_items(p.origin, 0);
    const double s = v->dtype == VIAL_DT_SET ? 1.0 : length_scale(arg_type(k, i));
    float*    xyz = md_alloc(k->alloc, (n ? n : 1) * 3 * sizeof(float));
    uint32_t* id = ids ? md_alloc(k->alloc, (n ? n : 1) * sizeof(uint32_t)) : NULL;
    const float* src = (const float*)p.data;
    size_t m = 0;
    for (size_t q = 0; q < n; ++q) {
        if (vial_val_is_missing(&p, 3 * q) || vial_val_is_missing(&p, 3 * q + 1) || vial_val_is_missing(&p, 3 * q + 2)) continue;
        for (int c = 0; c < 3; ++c) xyz[3 * m + c] = (float)(src[3 * q + c] * s);
        if (id) id[m] = ids[q];
        m += 1;
    }
    *out = (rdf_points_t){.xyz = xyz, .ids = id, .n = m};
    return true;
}

// How many points of a are also points of b: the pairs of a particle with itself, which are not counted
static uint64_t shared(vial_kctx_t* k, const rdf_points_t* a, const rdf_points_t* b) {
    if (!a->ids || !b->ids || !a->n || !b->n) return 0;
    uint32_t top = 0;
    for (size_t i = 0; i < b->n; ++i) if (b->ids[i] > top) top = b->ids[i];
    uint8_t* in_b = md_alloc(k->alloc, (size_t)top + 1);
    md_memset(in_b, 0, (size_t)top + 1);
    for (size_t i = 0; i < b->n; ++i) in_b[b->ids[i]] = 1;
    uint64_t n = 0;
    for (size_t i = 0; i < a->n; ++i) n += a->ids[i] <= top && in_b[a->ids[i]];
    return n;
}

// ### COUNTING ###

typedef struct rdf_count_t {
    const uint32_t* ida;
    const uint32_t* idb;
    float           lo2, hi2;
    double          lo, inv_width;
    uint32_t        bins;
    uint32_t*       counts;
} rdf_count_t;

static void count_pairs(const uint32_t* i_idx, const uint32_t* j_idx, const float* d2, size_t num_pairs, void* user) {
    rdf_count_t* c = user;
    for (size_t q = 0; q < num_pairs; ++q) {
        if (!(d2[q] >= c->lo2 && d2[q] < c->hi2)) continue;
        if (c->ida && c->idb && c->ida[i_idx[q]] == c->idb[j_idx[q]]) continue;
        const double b = (sqrt((double)d2[q]) - c->lo) * c->inv_width;
        if (!(b >= 0.0)) continue;
        c->counts[b < (double)c->bins ? (uint32_t)b : c->bins - 1] += 1;
    }
}

// ### THE ACCUMULATOR ###

static void free_frames(vial_rdf_t* rdf) {
    md_allocator_i* heap = md_get_heap_allocator();
    if (rdf->frame_counts) md_free(heap, rdf->frame_counts, (size_t)rdf->num_frames * rdf->bins * sizeof(uint32_t));
    if (rdf->frame_pairs) md_free(heap, rdf->frame_pairs, (size_t)rdf->num_frames * sizeof(double));
    rdf->frame_counts = NULL;
    rdf->frame_pairs = NULL;
}

// The sum, and what is kept per frame, set up for num_frames frames, nothing added; the counts per frame too when they
// fit in budget bytes, which keep gives what they take
static bool alloc_sums(vial_rdf_t* rdf, size_t budget, size_t* keep) {
    md_allocator_i* heap = md_get_heap_allocator();
    rdf->sum = md_alloc(heap, rdf->bins * sizeof(uint64_t));
    if (!rdf->sum) return false;
    md_memset(rdf->sum, 0, rdf->bins * sizeof(uint64_t));
    const size_t nf = rdf->num_frames > 0 ? (size_t)rdf->num_frames : 0;
    if (nf) {
        rdf->frame_pairs = md_alloc(heap, nf * sizeof(double));
        if (!rdf->frame_pairs) return false;
        for (size_t f = 0; f < nf; ++f) rdf->frame_pairs[f] = NAN;    // not added
        // Written only for the frames added, which frame_pairs tells apart, so not zeroed
        const size_t bytes = nf * rdf->bins * sizeof(uint32_t);
        if (keep && bytes / rdf->bins / sizeof(uint32_t) == nf && bytes <= budget) {
            rdf->frame_counts = md_alloc(heap, bytes);
            if (rdf->frame_counts) *keep = bytes;
        }
    }
    rdf->have_lock = md_mutex_init(&rdf->lock);
    return true;
}

bool vial_rdf_begin(vial_kctx_t* k, vial_rdf_t* rdf, const vial_args_t* a) {
    md_memset(rdf, 0, sizeof(*rdf));
    const int ci = vial_arg_index(a, "cutoff", 2);
    const vial_type_t* ct = arg_type(k, ci);
    double s1 = length_scale(ct), s2 = s1;
    if (ct && ct->dtype == VIAL_DT_TUPLE && ct->num_elems == 2) {
        s1 = length_scale(&ct->elems[0]);
        s2 = length_scale(&ct->elems[1]);
    }
    if (!vial_rdf_range(ci < 0 ? NULL : a->v[ci], s1, s2, &rdf->lo, &rdf->hi)) {
        return vial_rt_fail(k, "R0103", "cutoff= is a length more than 0, or lo..hi with hi above lo");
    }
    // The bins as the check worked them out: the extent of the type
    const vial_type_t* t = k->type;
    rdf->bins = t && t->ndim == 1 && t->dims[0].kind == VIAL_EXT_STATIC ? t->dims[0].n : 0;
    if (!rdf->bins) return vial_rt_fail(k, "R0103", "`rdf` has no bins");
    if (!vial_k_room(k, rdf->bins)) return false;
    rdf->num_frames = k->in->num_frames;
    vial_interp_t* in = k->in;
    size_t kept = 0;
    const bool ok = alloc_sums(rdf, in->frame_budget > in->frame_kept ? in->frame_budget - in->frame_kept : 0, &kept);
    in->frame_kept += kept;
    if (!ok) {
        const uint32_t bins = rdf->bins;
        vial_rdf_free(rdf);
        return vial_rt_fail(k, "R0105", "an rdf of %u bins does not fit in memory", bins);
    }
    return true;
}

bool vial_rdf_add(vial_kctx_t* k, vial_rdf_t* rdf, const vial_args_t* a) {
    const vial_frame_t* fr = k->in->frame;
    if (!fr || !fr->xyz) return vial_rt_fail(k, "R0103", "`rdf` adds frames");
    const int64_t f = fr->index;
    if (f < 0 || f >= rdf->num_frames) return vial_rt_fail(k, "R0103", "frame %lld is past the run's %lld", (long long)f, (long long)rdf->num_frames);
    if (!has_cell(&fr->cell)) {
        vial_rt_note(k, "R0101", "frame %lld has no unit cell; `rdf` leaves it out", (long long)f);
        return true;
    }
    double width;
    const double vol = cell_volume_width(&fr->cell, &width);
    if (!(vol > 0.0)) {
        vial_rt_note(k, "R0101", "frame %lld has a cell of no volume; `rdf` leaves it out", (long long)f);
        return true;
    }
    if (2.0 * rdf->hi > width) {
        vial_rt_note(k, "R0101", "frame %lld: the cell is %.4g Å across at its narrowest, less than twice the cutoff, %.4g Å; `rdf` leaves it out",
            (long long)f, width, rdf->hi);
        return true;
    }
    rdf_points_t pa, pb;
    if (!points_of(k, a, "a", 0, &pa) || !points_of(k, a, "b", 1, &pb)) return false;
    uint32_t* counts = md_alloc(k->alloc, rdf->bins * sizeof(uint32_t));
    md_memset(counts, 0, rdf->bins * sizeof(uint32_t));
    if (pa.n && pb.n) {
        md_spatial_acc_t acc = {.alloc = k->alloc};
        md_spatial_acc_init(&acc, &(md_spatial_acc_desc_t){.xyz = pb.xyz, .count = pb.n, .cutoff = rdf->hi, .unitcell = &fr->cell});
        rdf_count_t c = {
            .ida = pa.ids, .idb = pb.ids,
            .lo2 = (float)(rdf->lo * rdf->lo), .hi2 = (float)(rdf->hi * rdf->hi),
            .lo = rdf->lo, .inv_width = rdf->bins / (rdf->hi - rdf->lo),
            .bins = rdf->bins, .counts = counts,
        };
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, pa.xyz, 0, NULL, pa.n, rdf->hi, count_pairs, &c, 0);
        md_spatial_acc_free(&acc);
    }
    // The pairs there are, spread evenly over the cell: per Å³
    const double pairs = ((double)pa.n * (double)pb.n - (double)shared(k, &pa, &pb)) / vol;
    // Each frame its own entries, so no lock for those; the sum under it
    if (rdf->frame_counts) md_memcpy(rdf->frame_counts + (size_t)f * rdf->bins, counts, rdf->bins * sizeof(uint32_t));
    rdf->frame_pairs[f] = pairs;
    if (rdf->have_lock) md_mutex_lock(&rdf->lock);
    for (uint32_t b = 0; b < rdf->bins; ++b) rdf->sum[b] += counts[b];
    if (rdf->have_lock) md_mutex_unlock(&rdf->lock);
    return true;
}

// The rdf from counts of pairs and P, the pairs per Å³ summed over the same frames; NaN in every bin when P is 0
static void normalize(float* out, const uint64_t* counts, double P, const vial_rdf_t* rdf) {
    const double w = (rdf->hi - rdf->lo) / rdf->bins;
    for (uint32_t b = 0; b < rdf->bins; ++b) {
        const double r0 = rdf->lo + b * w, r1 = rdf->lo + (b + 1) * w;
        const double shell = 4.0 / 3.0 * MD_PI * (r1 * r1 * r1 - r0 * r0 * r0);
        out[b] = P > 0.0 ? (float)((double)counts[b] / (P * shell)) : NAN;
    }
}

// The pairs per Å³ over the frames added where mask is set (all when NULL), added in frame order
static double pairs_over(const vial_rdf_t* rdf, const uint8_t* mask) {
    double P = 0.0;
    for (int64_t f = 0; f < rdf->num_frames; ++f) {
        if (mask && !mask[f]) continue;
        if (!isnan(rdf->frame_pairs[f])) P += rdf->frame_pairs[f];
    }
    return P;
}

bool vial_rdf_end(vial_kctx_t* k, vial_rdf_t* rdf, vial_val_t* out) {
    vial_val_t o;
    if (!vial_k_make(k, VIAL_DT_FLOAT, 1, &rdf->bins, &o)) return false;
    float* g = (float*)o.data;
    normalize(g, rdf->sum, rdf->frame_pairs ? pairs_over(rdf, NULL) : 0.0, rdf);
    for (uint32_t b = 0; b < rdf->bins; ++b) if (isnan(g[b])) vial_val_set_missing(&o, b, k->alloc);
    *out = o;
    return true;
}

bool vial_rdf_over(const vial_rdf_t* rdf, const uint8_t* mask, float* out) {
    if (!rdf || !rdf->frame_counts || !rdf->frame_pairs || !mask) return false;
    uint64_t* counts = md_alloc(md_get_heap_allocator(), rdf->bins * sizeof(uint64_t));
    if (!counts) return false;
    md_memset(counts, 0, rdf->bins * sizeof(uint64_t));
    for (int64_t f = 0; f < rdf->num_frames; ++f) {
        if (!mask[f] || isnan(rdf->frame_pairs[f])) continue;
        const uint32_t* c = rdf->frame_counts + (size_t)f * rdf->bins;
        for (uint32_t b = 0; b < rdf->bins; ++b) counts[b] += c[b];
    }
    normalize(out, counts, pairs_over(rdf, mask), rdf);
    md_free(md_get_heap_allocator(), counts, rdf->bins * sizeof(uint64_t));
    return true;
}

bool vial_rdf_like(vial_rdf_t* out, const vial_rdf_t* rdf) {
    md_memset(out, 0, sizeof(*out));
    out->bins = rdf->bins;
    out->lo = rdf->lo;
    out->hi = rdf->hi;
    out->num_frames = rdf->num_frames;
    if (alloc_sums(out, 0, NULL)) return true;
    vial_rdf_free(out);
    return false;
}

void vial_rdf_free(vial_rdf_t* rdf) {
    if (!rdf) return;
    if (rdf->sum) md_free(md_get_heap_allocator(), rdf->sum, rdf->bins * sizeof(uint64_t));
    free_frames(rdf);
    if (rdf->have_lock) md_mutex_destroy(&rdf->lock);
    md_memset(rdf, 0, sizeof(*rdf));
}

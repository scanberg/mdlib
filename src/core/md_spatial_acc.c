#include "md_spatial_acc.h"

#include <core/md_coord_stream.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_intrinsics.h>
#include <md_util.h>
#include <md_types.h>

#include <float.h>
#include <math.h>

#define SPATIAL_ACC_BUFLEN 1024
#define SPATIAL_ACC_MAX_NEIGHBOR_CELLS 5

typedef md_128i ivec4_t;

static inline ivec4_t ivec4_set(int x, int y, int z, int w) {
    return simde_mm_set_epi32(w, z, y, x);
}

static inline ivec4_t ivec4_load(const int* v) {
    return simde_mm_loadu_si128((const md_128i*)v);
}

static inline void ivec4_store(int* v, ivec4_t a) {
    simde_mm_storeu_si128((md_128i*)v, a);
}

static inline ivec4_t ivec4_set1(int v) {
    return simde_mm_set1_epi32(v);
}

static inline ivec4_t ivec4_min(ivec4_t a, ivec4_t b) {
    return simde_mm_min_epi32(a, b);
}

static inline ivec4_t ivec4_max(ivec4_t a, ivec4_t b) {
    return simde_mm_max_epi32(a, b);
}

static inline ivec4_t ivec4_clamp(ivec4_t v, ivec4_t min, ivec4_t max) {
    return ivec4_max(ivec4_min(v, max), min);
}

static inline ivec4_t ivec4_from_vec4(vec4_t v) {
    return simde_mm_cvtps_epi32(v.m128);
}

static inline vec4_t vec4_from_ivec4(ivec4_t v) {
    vec4_t r;
    r.m128 = md_mm_cvtepi32_ps(v);
    return r;
}

static inline ivec4_t ivec4_add(ivec4_t a, ivec4_t b) {
    return simde_mm_add_epi32(a, b);
}

static inline ivec4_t ivec4_sub(ivec4_t a, ivec4_t b) {
    return simde_mm_sub_epi32(a, b);
}

static inline ivec4_t ivec4_cmpgt(ivec4_t a, ivec4_t b) {
    return simde_mm_cmpgt_epi32(a, b);
}

static inline ivec4_t ivec4_cmplt(ivec4_t a, ivec4_t b) {
    return simde_mm_cmplt_epi32(a, b);
}

static inline ivec4_t ivec4_and(ivec4_t a, ivec4_t b) {
    return simde_mm_and_si128(a, b);
}

static inline ivec4_t ivec4_andnot(ivec4_t a, ivec4_t b) {
    return simde_mm_andnot_si128(b, a);
}

static inline ivec4_t ivec4_or(ivec4_t a, ivec4_t b) {
    return simde_mm_or_si128(a, b);
}

static inline bool ivec4_any(ivec4_t v) {
#if defined(__aarch64__)
    // Use NEON reduction: true if any 32-bit lane is non-zero.
    simde__m128i_private vp = simde__m128i_to_private(v);
    return vmaxvq_u32(vp.neon_u32) != 0u;
#else
    // Efficient on x86: ptest(a,a) via SIMDe
    return !simde_mm_testz_si128(v, v);
#endif
}

typedef struct {
    float x, y, z;
    uint32_t idx;
} elem_t;

static inline vec4_t md_coord_stream_load_vec4(const md_coord_stream_t* s, size_t i) {
    const size_t src = s->idx ? (size_t)s->idx[i] : i;

    if (s->layout == MD_COORD_STREAM_LAYOUT_SOA) {
        return vec4_set(s->soa.x[src], s->soa.y[src], s->soa.z[src], 0.0f);
    } else {
        const char* p = (const char*)s->aos.base + src * s->aos.stride;
        vec4_t v = {0};
        MEMCPY(&v, p, sizeof(float) * 3);
        return v;
    }
}

static inline int32_t md_coord_stream_load_idx(const md_coord_stream_t* s, size_t i) {
    return s->idx ? s->idx[i] : (int32_t)i;
}

void md_spatial_acc_free(md_spatial_acc_t* acc) {
    ASSERT(acc);
    if (acc->alloc) {
        if (acc->elem_x)   md_array_free(acc->elem_x,   acc->alloc);
        if (acc->elem_y)   md_array_free(acc->elem_y,   acc->alloc);
        if (acc->elem_z)   md_array_free(acc->elem_z,   acc->alloc);
        if (acc->elem_idx) md_array_free(acc->elem_idx, acc->alloc);
        if (acc->cell_off) md_array_free(acc->cell_off, acc->alloc);
        if (acc->top_mask)       md_array_free(acc->top_mask,       acc->alloc);
        if (acc->top_base)       md_array_free(acc->top_base,       acc->alloc);
        for (int t = 0; t < MD_SPATIAL_ACC_MAX_TIERS; ++t) {
            if (acc->tier_mask[t]) md_array_free(acc->tier_mask[t], acc->alloc);
            if (acc->tier_base[t]) md_array_free(acc->tier_base[t], acc->alloc);
        }
    }
    MEMSET(acc, 0, sizeof(md_spatial_acc_t));
}

static void md_spatial_acc_reset(md_spatial_acc_t* acc) {
    ASSERT(acc);
    md_array_shrink(acc->elem_x, 0);
	md_array_shrink(acc->elem_y, 0);
	md_array_shrink(acc->elem_z, 0);
	md_array_shrink(acc->elem_idx, 0);
	md_array_shrink(acc->cell_off, 0);
	md_array_shrink(acc->top_mask, 0);
	md_array_shrink(acc->top_base, 0);
    for (int t = 0; t < MD_SPATIAL_ACC_MAX_TIERS; ++t) {
        md_array_shrink(acc->tier_mask[t], 0);
        md_array_shrink(acc->tier_base[t], 0);
    }
    acc->num_tiers = 0;
    MEMSET(acc->top_dim, 0, sizeof(acc->top_dim));

    acc->num_cells = 0;
    acc->num_elems = 0;

    acc->G00 = acc->G11 = acc->G22 = 0.0f;
    acc->H01 = acc->H02 = acc->H12 = 0.0f;

    MEMSET(acc->cell_dim,  0, sizeof(acc->cell_dim));

	acc->flags = 0;
	MEMSET(acc->origin, 0, sizeof(acc->origin));
}

// ### CELL INDEX ###
//
// The elements are binned into cells and sorted by cell. The cells themselves are not a dense grid, which would store
// an offset for every cell of the frame whether occupied or not: its size would follow the volume of the frame over
// the cell volume, whatever the number of elements, and a large sparse system (fibrils in a box of micrometres) or a
// small subset of a large system would cost gigabytes, mostly for empty cells.
//
// Only the occupied cells are stored, in a hierarchy of 4x4x4 groups:
//   tier 0           the occupied cells: cell_off[k] .. cell_off[k+1] are the elements of occupied cell k
//   tier 1 .. L-1    the occupied nodes, each with a 64 bit mask of its occupied children and the index of the first
//   tier L (top)     a dense grid of nodes over the frame, 4^L cells along each axis per node
// The child at local position (x, y, z) in 0..3 is bit x | y << 2 | z << 4 of its parent's mask, and its index is
// base + popcount(mask & ((1 << bit) - 1)). Finding a cell is L steps of a load and a popcount, no search and no hash.
// The cells (and the elements) are ordered top node first, then depth first through the tiers in bit order, so the
// cells of every node are contiguous, x running fastest within a node.
//
// L is the smallest number of tiers (at least 2) for which the top grid has no more nodes than max(N, 2^16), so the
// only part which grows with the volume of the frame is bounded by the number of elements. Everything else is
// proportional to the number of occupied cells.

#define SPATIAL_ACC_MAX_CELLS_PER_DIM (1u << 20)
#define SPATIAL_ACC_MIN_TOP_LIMIT     (1u << 16)

// Local position of the tier t node containing cell c, within its parent (the tier t+1 node): the bit of its parent's mask
static inline uint32_t local_bit(const uint32_t c[3], uint32_t t) {
    const uint32_t s = 2 * t;
    return ((c[0] >> s) & 3) | (((c[1] >> s) & 3) << 2) | (((c[2] >> s) & 3) << 4);
}

static inline size_t top_index(const md_spatial_acc_t* acc, const uint32_t c[3]) {
    const uint32_t s = 2 * acc->num_tiers;
    return ((size_t)(c[2] >> s) * acc->top_dim[1] + (size_t)(c[1] >> s)) * acc->top_dim[0] + (size_t)(c[0] >> s);
}

// The node of tier T (1 .. num_tiers) containing cell coordinates c (within the grid): its mask and first child, from
// the top down. A mask of 0 if there is no such node.
static inline void node_lookup(const md_spatial_acc_t* acc, const uint32_t c[3], uint32_t T, uint64_t* out_mask, uint32_t* out_base) {
    const size_t t = top_index(acc, c);
    uint64_t mask = acc->top_mask[t];
    uint32_t base = acc->top_base[t];
    for (uint32_t l = acc->num_tiers; l > T; --l) {
        const uint32_t bit = local_bit(c, l - 1);
        if (!((mask >> bit) & 1)) {
            *out_mask = 0;
            *out_base = 0;
            return;
        }
        const uint32_t idx = base + (uint32_t)popcnt64(mask & ((1ULL << bit) - 1));
        mask = acc->tier_mask[l - 1][idx];
        base = acc->tier_base[l - 1][idx];
    }
    *out_mask = mask;
    *out_base = base;
}

// Index of the occupied cell at cell coordinates c (within the grid), or UINT32_MAX if it is empty
static inline uint32_t cell_lookup(const md_spatial_acc_t* acc, const uint32_t c[3]) {
    uint64_t mask;
    uint32_t base;
    node_lookup(acc, c, 1, &mask, &base);
    const uint32_t bit = local_bit(c, 0);
    if (!((mask >> bit) & 1)) return UINT32_MAX;
    return base + (uint32_t)popcnt64(mask & ((1ULL << bit) - 1));
}

// Cell coordinates of a fractional coordinate, exactly as the elements were binned
static inline void bin_cell(uint32_t out[3], float sx, float sy, float sz, vec4_t fcell_dim, ivec4_t icell_max) {
    const vec4_t s = vec4_set(sx, sy, sz, 0);
    ivec4_t ic = ivec4_from_vec4(vec4_floor(vec4_mul(s, fcell_dim)));
    ic = ivec4_clamp(ic, ivec4_set1(0), icell_max);
    uint32_t tmp[4];
    md_mm_storeu_epi32(tmp, ic);
    out[0] = tmp[0];
    out[1] = tmp[1];
    out[2] = tmp[2];
}

// Builds the elements and the cell index. The frame (cell_dim, origin, I, flags) is already stored in acc.
static void build_cells(md_spatial_acc_t* acc, const md_coord_stream_t* coords, md_spatial_acc_flags_t in_flags) {
    const size_t N = coords->count;
    const uint32_t* cell_dim = acc->cell_dim;

    // Number of tiers: the top grid has at most max(N, 2^16) nodes
    const size_t top_limit = MAX(N, (size_t)SPATIAL_ACC_MIN_TOP_LIMIT);
    uint32_t L = 2;
    uint32_t top_dim[3];
    size_t num_top;
    for (;;) {
        const uint32_t s = 2 * L;
        for (int i = 0; i < 3; ++i) {
            top_dim[i] = (uint32_t)(((uint64_t)cell_dim[i] + (1ULL << s) - 1) >> s);
        }
        num_top = (size_t)top_dim[0] * top_dim[1] * top_dim[2];
        if (num_top <= top_limit || L == MD_SPATIAL_ACC_MAX_TIERS) break;
        ++L;
    }
    acc->num_tiers = L;
    MEMCPY(acc->top_dim, top_dim, sizeof(top_dim));

    md_temp_scope_t temp_scope = md_temp_begin_avoid(acc->alloc);
    elem_t*   buf[2] = {
        (elem_t*)md_temp_alloc(temp_scope, MAX(N, 1) * sizeof(elem_t)),
        (elem_t*)md_temp_alloc(temp_scope, MAX(N, 1) * sizeof(elem_t)),
    };
    uint32_t* digit = (uint32_t*)md_temp_alloc(temp_scope, MAX(N, 1) * sizeof(uint32_t));
    const size_t num_hist = MAX(num_top, 4096) + 1;
    uint32_t* hist = (uint32_t*)md_temp_alloc(temp_scope, num_hist * sizeof(uint32_t));

    const vec4_t  fcell_dim = vec4_set((float)cell_dim[0], (float)cell_dim[1], (float)cell_dim[2], 0);
    const ivec4_t icell_max = ivec4_set(cell_dim[0] - 1, cell_dim[1] - 1, cell_dim[2] - 1, 0);

    float val;
    MEMSET(&val, 0xFF, sizeof(val));
    const uint32_t flags = acc->flags;
    const vec4_t pbc_mask = vec4_set((flags & MD_UNITCELL_PBC_X) ? val : 0, (flags & MD_UNITCELL_PBC_Y) ? val : 0, (flags & MD_UNITCELL_PBC_Z) ? val : 0, 0);
    const vec4_t origin = vec4_set(acc->origin[0], acc->origin[1], acc->origin[2], 0);
    const vec4_t vI[3] = {
        vec4_set(acc->I[0][0], acc->I[0][1], acc->I[0][2], 0),
        vec4_set(acc->I[1][0], acc->I[1][1], acc->I[1][2], 0),
        vec4_set(acc->I[2][0], acc->I[2][1], acc->I[2][2], 0),
    };

    // 1) Fractional coordinates, periodic axes wrapped into [0,1). The record keeps the position in the stream.
    for (size_t i = 0; i < N; ++i) {
        const vec4_t r = md_coord_stream_load_vec4(coords, i);
        vec4_t s = vec4_linear_combine_3(vec4_sub(r, origin), vI);
        s = vec4_blend(s, vec4_fract(s), pbc_mask);
        buf[0][i] = (elem_t){ s.x, s.y, s.z, (uint32_t)i };
    }

    // 2) Stable counting sorts, least significant first: the local positions two tiers at a time, then the top node
    int src = 0;
    const uint32_t num_local_passes = (L + 1) / 2;
    for (uint32_t pass = 0; pass <= num_local_passes; ++pass) {
        const bool top_pass = pass == num_local_passes;
        const size_t num_buckets = top_pass ? num_top : 4096;
        MEMSET(hist, 0, (num_buckets + 1) * sizeof(uint32_t));

        for (size_t i = 0; i < N; ++i) {
            uint32_t c[3];
            bin_cell(c, buf[src][i].x, buf[src][i].y, buf[src][i].z, fcell_dim, icell_max);
            uint32_t d;
            if (top_pass) {
                d = (uint32_t)top_index(acc, c);
            } else {
                const uint32_t t0 = 2 * pass;
                d = local_bit(c, t0);
                if (t0 + 1 < L) d |= local_bit(c, t0 + 1) << 6;
            }
            digit[i] = d;
            hist[d + 1] += 1;
        }
        for (size_t b = 0; b < num_buckets; ++b) {
            hist[b + 1] += hist[b];
        }

        if (!top_pass) {
            elem_t* dst = buf[src ^ 1];
            for (size_t i = 0; i < N; ++i) {
                dst[hist[digit[i]]++] = buf[src][i];
            }
            src ^= 1;
        } else {
            // The last pass scatters into the elements
            const size_t alloc_len = ALIGN_TO(N + 8, 16);
            md_array_resize(acc->elem_x,   alloc_len, acc->alloc);
            md_array_resize(acc->elem_y,   alloc_len, acc->alloc);
            md_array_resize(acc->elem_z,   alloc_len, acc->alloc);
            md_array_resize(acc->elem_idx, alloc_len, acc->alloc);
            // The padding is read (masked) by the vectorized loops
            MEMSET(acc->elem_x + N, 0, (alloc_len - N) * sizeof(float));
            MEMSET(acc->elem_y + N, 0, (alloc_len - N) * sizeof(float));
            MEMSET(acc->elem_z + N, 0, (alloc_len - N) * sizeof(float));
            MEMSET(acc->elem_idx + N, 0, (alloc_len - N) * sizeof(uint32_t));
            for (size_t i = 0; i < N; ++i) {
                const elem_t e = buf[src][i];
                const uint32_t dst = hist[digit[i]]++;
                acc->elem_x[dst] = e.x;
                acc->elem_y[dst] = e.y;
                acc->elem_z[dst] = e.z;
                const uint32_t stream_idx = (uint32_t)md_coord_stream_load_idx(coords, e.idx);
                acc->elem_idx[dst] = (in_flags & MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX) ? stream_idx : e.idx;
            }
        }
    }
    acc->num_elems = N;

    // 3) Occupied cells and the tiers, in one pass over the sorted elements. A cell which differs from the previous
    //    one in the node of tier t (or above) starts new nodes in the tiers t-1 .. 1, and a new top node when it
    //    differs from it in the top node. Each node records the index of its first child as its base when it starts.
    md_array_resize(acc->top_mask, num_top, acc->alloc);
    md_array_resize(acc->top_base, num_top, acc->alloc);
    MEMSET(acc->top_mask, 0, num_top * sizeof(uint64_t));
    MEMSET(acc->top_base, 0, num_top * sizeof(uint32_t));
    md_array_shrink(acc->cell_off, 0);

    uint32_t prev[3] = {0};
    size_t   top_idx = 0;
    uint32_t num_occ = 0;
    for (size_t k = 0; k < N; ++k) {
        uint32_t c[3];
        bin_cell(c, acc->elem_x[k], acc->elem_y[k], acc->elem_z[k], fcell_dim, icell_max);
        if (k > 0 && c[0] == prev[0] && c[1] == prev[1] && c[2] == prev[2]) continue;

        // The highest tier whose node changes (L for the top node, 0 if only the cell changes)
        uint32_t hi = L;
        if (k > 0) {
            hi = 0;
            for (uint32_t t = L; t >= 1; --t) {
                const uint32_t s = 2 * t;
                if ((c[0] >> s) != (prev[0] >> s) || (c[1] >> s) != (prev[1] >> s) || (c[2] >> s) != (prev[2] >> s)) {
                    hi = t;
                    break;
                }
            }
        }
        if (hi == L) {
            top_idx = top_index(acc, c);
            acc->top_base[top_idx] = (uint32_t)md_array_size(acc->tier_mask[L - 1]);
        }
        for (uint32_t t = MIN(hi, L - 1); t >= 1; --t) {
            const uint32_t first_child = (t == 1) ? num_occ : (uint32_t)md_array_size(acc->tier_mask[t - 1]);
            md_array_push(acc->tier_mask[t], 0, acc->alloc);
            md_array_push(acc->tier_base[t], first_child, acc->alloc);
        }

        acc->top_mask[top_idx] |= 1ULL << local_bit(c, L - 1);
        for (uint32_t t = 1; t < L; ++t) {
            acc->tier_mask[t][md_array_size(acc->tier_mask[t]) - 1] |= 1ULL << local_bit(c, t - 1);
        }

        md_array_push(acc->cell_off, (uint32_t)k, acc->alloc);
        num_occ += 1;
        MEMCPY(prev, c, sizeof(prev));
    }
    md_array_push(acc->cell_off, (uint32_t)N, acc->alloc);
    acc->num_cells = num_occ;

    md_temp_end(temp_scope);
}

// ### CELL EXTENT FROM OCCUPANCY ###
//
// The pair loops go cell by cell, testing the points of a cell against those of its neighbours eight at a time. When
// the cells hold a point or two, most of that work is overhead: the cell lookups, and SIMD lanes with nothing in them.
// For a cutoff of 10 A, 4 million coarse grained beads in a 2 x 2 x 8 um box sit about one per 10 A cell, and the
// query runs 4x faster with 40 A cells. Points at atomistic density fill 10 A cells with ~100 and want the cutoff.
//
// So the cell extent is the smallest of a ladder of multiples of the cutoff (steps of 2^(1/4)) at which a point finds, on average,
// SPATIAL_ACC_TARGET_OCCUPANCY points in its own cell. That is weighted by points rather than by cells, so it follows
// where the points are, not how much empty space there is between them. It is estimated from a random subset of the
// points (each kept with probability f): a point which shares its cell with k - 1 others of the subset sees about
// 1 + (k - 1) / f points of the whole set there.

#ifndef SPATIAL_ACC_TARGET_OCCUPANCY
#define SPATIAL_ACC_TARGET_OCCUPANCY 10.0
#endif
#define SPATIAL_ACC_MAX_CELL_SCALE   4.0
#define SPATIAL_ACC_DEFAULT_CUTOFF   6.0     // For a description without a cutoff
#define SPATIAL_ACC_OCCUPANCY_SAMPLES 16384

static inline uint64_t occupancy_hash(uint64_t x) {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

static void sort_u64(uint64_t* keys, uint64_t* tmp, size_t n) {
    // LSD radix, 8 bits at a time, passes where every key has the same digit are skipped
    uint32_t hist[256];
    for (int shift = 0; shift < 64; shift += 8) {
        MEMSET(hist, 0, sizeof(hist));
        for (size_t i = 0; i < n; ++i) hist[(keys[i] >> shift) & 0xFF] += 1;
        if (n == 0 || hist[(keys[0] >> shift) & 0xFF] == n) continue;
        uint32_t sum = 0;
        for (int b = 0; b < 256; ++b) {
            const uint32_t c = hist[b];
            hist[b] = sum;
            sum += c;
        }
        for (size_t i = 0; i < n; ++i) tmp[hist[(keys[i] >> shift) & 0xFF]++] = keys[i];
        MEMCPY(keys, tmp, n * sizeof(uint64_t));
    }
}

static double cell_ext_from_occupancy(const md_coord_stream_t* coords, double cutoff, const md_allocator_i* avoid) {
    const size_t N = coords->count;
    if (N == 0 || !(cutoff > 0.0)) return cutoff;

    md_temp_scope_t temp = md_temp_begin_avoid(avoid);
    const size_t cap = MIN(N, (size_t)SPATIAL_ACC_OCCUPANCY_SAMPLES * 2);
    vec4_t*   pos  = md_temp_alloc_array(temp, vec4_t, cap);
    uint64_t* key  = md_temp_alloc_array(temp, uint64_t, cap);
    uint64_t* tmp  = md_temp_alloc_array(temp, uint64_t, cap);

    // The subset: every point when there are few, otherwise each with probability f, decided by a hash of its position
    // in the stream. Positions in a stream follow molecules, so taking every n-th would not be a random subset.
    const double f = N <= SPATIAL_ACC_OCCUPANCY_SAMPLES ? 1.0 : (double)SPATIAL_ACC_OCCUPANCY_SAMPLES / (double)N;
    const uint64_t threshold = f >= 1.0 ? UINT64_MAX : (uint64_t)(f * 18446744073709551616.0);
    size_t m = 0;
    vec4_t lo = vec4_set1(FLT_MAX);
    for (size_t i = 0; i < N && m < cap; ++i) {
        if (f < 1.0 && occupancy_hash(i) >= threshold) continue;
        const vec4_t p = md_coord_stream_load_vec4(coords, i);
        pos[m++] = p;
        lo = vec4_min(lo, p);
    }

    double result = cutoff * SPATIAL_ACC_MAX_CELL_SCALE;
    // Steps of 2^(1/4)
    for (double scale = 1.0; scale <= SPATIAL_ACC_MAX_CELL_SCALE * 1.0001; scale *= 1.189207115002721) {
        const double ext = cutoff * scale;
        const float inv = (float)(1.0 / ext);
        for (size_t i = 0; i < m; ++i) {
            const vec4_t c = vec4_mul1(vec4_sub(pos[i], lo), inv);
            const uint64_t cx = (uint64_t)MIN(c.x, 2097151.0f);
            const uint64_t cy = (uint64_t)MIN(c.y, 2097151.0f);
            const uint64_t cz = (uint64_t)MIN(c.z, 2097151.0f);
            key[i] = (cx << 42) | (cy << 21) | cz;
        }
        sort_u64(key, tmp, m);
        double sum = 0.0;
        for (size_t i = 0; i < m;) {
            size_t j = i + 1;
            while (j < m && key[j] == key[i]) ++j;
            const double k = (double)(j - i);
            sum += k * (1.0 + (k - 1.0) / f);
            i = j;
        }
        if (sum / (double)m >= SPATIAL_ACC_TARGET_OCCUPANCY) {
            result = ext;
            break;
        }
    }

    md_temp_end(temp);
    return result;
}

// The frame of the structure: grid dimensions, metric, basis and origin
// A and I are not const: C does not convert double (*)[3] to const double (*)[3] implicitly
static void store_frame(md_spatial_acc_t* acc, double A[3][3], double I[3][3], vec4_t origin, uint32_t flags,
                        const uint32_t cell_dim[3], const float inv_cell_ext[3],
                        double G00, double G11, double G22, double H01, double H02, double H12) {
    MEMCPY(acc->cell_dim,  cell_dim,  sizeof(acc->cell_dim));
    MEMCPY(acc->inv_cell_ext, inv_cell_ext, sizeof(acc->inv_cell_ext));

    acc->G00 = (float)G00;
    acc->G11 = (float)G11;
    acc->G22 = (float)G22;
    acc->H01 = (float)H01;
    acc->H02 = (float)H02;
    acc->H12 = (float)H12;

	acc->A[0][0] = (float)A[0][0];
	acc->A[0][1] = (float)A[0][1];
	acc->A[0][2] = (float)A[0][2];

	acc->A[1][0] = (float)A[1][0];
	acc->A[1][1] = (float)A[1][1];
	acc->A[1][2] = (float)A[1][2];

	acc->A[2][0] = (float)A[2][0];
	acc->A[2][1] = (float)A[2][1];
	acc->A[2][2] = (float)A[2][2];

    acc->I[0][0] = (float)I[0][0];
    acc->I[0][1] = (float)I[0][1];
    acc->I[0][2] = (float)I[0][2];

    acc->I[1][0] = (float)I[1][0];
    acc->I[1][1] = (float)I[1][1];
    acc->I[1][2] = (float)I[1][2];

    acc->I[2][0] = (float)I[2][0];
    acc->I[2][1] = (float)I[2][1];
    acc->I[2][2] = (float)I[2][2];

    // Persist origin offset used to construct fractional frame
    acc->origin[0] = origin.x;
    acc->origin[1] = origin.y;
    acc->origin[2] = origin.z;

    acc->flags = flags;
}

static void spatial_acc_init_internal(md_spatial_acc_t* acc, const md_coord_stream_t* coords, double in_cutoff, const md_unitcell_t* in_unitcell, md_spatial_acc_flags_t in_flags) {
    ASSERT(acc);
    ASSERT(coords);

    if (!acc->alloc) {
        MD_LOG_ERROR("Must have allocator set within spatial acc");
        return;
    }

    if (in_flags & MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX) {
        if (!coords->idx) {
            MD_LOG_ERROR("Flag MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX is set, but coordinate stream index is not supplied");
            return;
        }
    }

	// Reset acc, but to not free memory
    md_spatial_acc_reset(acc);

    if (coords->count == 0) {
		// Not a real error, but nothing to build. Leave acc in a valid empty state.
        return;
    }

    // The cells: at least the cutoff, so a query up to it reaches one cell around, larger where the points are sparse.
    // Never below 3 A, which would only multiply the cells to visit.
    const double cutoff   = in_cutoff > 0.0 ? in_cutoff : SPATIAL_ACC_DEFAULT_CUTOFF;
    const double CELL_EXT = MAX(cell_ext_from_occupancy(coords, cutoff, acc->alloc), 3.0);

    double A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double I[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    uint32_t flags = 0;

    if (in_unitcell) {
        md_unitcell_A_extract_double(A, in_unitcell);
        md_unitcell_I_extract_double(I, in_unitcell);
        flags = md_unitcell_flags(in_unitcell);
    }

    vec4_t origin = vec4_zero();

    if ((flags & MD_UNITCELL_PBC_ALL) != MD_UNITCELL_PBC_ALL) {
        ASSERT((flags & MD_UNITCELL_TRICLINIC) == 0);
        // Unit cell either missing or not periodic along one or more axis
        // Seed from the first point, not from the origin. A system which sits far from the origin would otherwise
        // get a grid spanning all the way back to it, and the cell count grows with that distance: two points at
        // 6700 A with an 8 A cell extent produce 840 cells per axis instead of 2.
        vec4_t aabb_min = md_coord_stream_load_vec4(coords, 0);
        vec4_t aabb_max = aabb_min;
        for (size_t i = 1; i < coords->count; i++) {
            vec4_t v = md_coord_stream_load_vec4(coords, i);
            aabb_min = vec4_min(aabb_min, v);
            aabb_max = vec4_max(aabb_max, v);
        }

        // Construct A and I from aabb extent
        vec4_t aabb_ext = vec4_sub(aabb_max, aabb_min);

        // Round up to nearest N * CELL_EXT
        aabb_ext = vec4_mul1(vec4_ceil(vec4_div1(aabb_ext, (float)CELL_EXT)), (float)CELL_EXT);

        // Set min as center - half extent
        vec4_t aabb_center = vec4_mul1(vec4_add(aabb_min, aabb_max), 0.5f);
        aabb_min = vec4_sub(aabb_center, vec4_mul1(aabb_ext, 0.5f));
        aabb_max = vec4_add(aabb_center, vec4_mul1(aabb_ext, 0.5f));

        if ((flags & MD_UNITCELL_PBC_X) == 0) {
            origin.x = aabb_min.x;
            if (aabb_ext.x > 0.0f) {
                A[0][0] = aabb_ext.x;
                I[0][0] = 1.0 / aabb_ext.x;
            }
        }
        if ((flags & MD_UNITCELL_PBC_Y) == 0) {
            origin.y = aabb_min.y;
            if (aabb_ext.y > 0.0f) {
                A[1][1] = aabb_ext.y;
                I[1][1] = 1.0 / aabb_ext.y;
            }
        }
        if ((flags & MD_UNITCELL_PBC_Z) == 0) {
            origin.z = aabb_min.z;
            if (aabb_ext.z > 0.0f) {
                A[2][2] = aabb_ext.z;
                I[2][2] = 1.0 / aabb_ext.z;
            }
        }
    }

    // Basis vectors are the COLUMNS of A (column-major storage: A[col][row])
    double a[3] = { A[0][0], A[0][1], A[0][2] };
    double b[3] = { A[1][0], A[1][1], A[1][2] };
    double c[3] = { A[2][0], A[2][1], A[2][2] };

    // Metric G = A^T A
    double G00 = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); // dot(a, a)
    double G11 = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]); // dot(b, b)
    double G22 = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]); // dot(c, c)
    double G01 = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]); // dot(a, b)
    double G02 = (a[0] * c[0] + a[1] * c[1] + a[2] * c[2]); // dot(a, c)
    double G12 = (b[0] * c[0] + b[1] * c[1] + b[2] * c[2]); // dot(b, c)

    double H01 = 0.0, H02 = 0.0, H12 = 0.0;

    double norm_a = sqrt(G00);
    double norm_b = sqrt(G11);
    double norm_c = sqrt(G22);

    float inv_cell_ext[3] = {
        (float)(norm_a > 0.0 ? 1.0 / norm_a : 0.0),
        (float)(norm_b > 0.0 ? 1.0 / norm_b : 0.0),
        (float)(norm_c > 0.0 ? 1.0 / norm_c : 0.0),
    };

    if (flags & MD_UNITCELL_TRICLINIC) {
        H01 = 2.0 * G01;
        H02 = 2.0 * G02;
        H12 = 2.0 * G12;

        double det = G00 * (G11 * G22 - G12 * G12) - G01 * (G01 * G22 - G12 * G02) + G02 * (G01 * G12 - G11 * G02);

        if (det < DBL_EPSILON) {
            MD_LOG_ERROR("Degenerate unit cell / A matrix provided to spatial acc");
            md_spatial_acc_reset(acc);
            return;
        }

        double GI00 = (G11 * G22 - G12 * G12) / det;
        double GI11 = (G00 * G22 - G02 * G02) / det;
        double GI22 = (G00 * G11 - G01 * G01) / det;

        inv_cell_ext[0] = (float)sqrt(GI00);
        inv_cell_ext[1] = (float)sqrt(GI11);
        inv_cell_ext[2] = (float)sqrt(GI22);
    }

    // Estimate cell_dim by measuring the extents of the box vectors (norms of columns of A)
    // This is only a heuristic for bin counts; the grid is still in fractional space.
    uint32_t cell_dim[3] = {
        (uint32_t)CLAMP(floor(norm_a / CELL_EXT), 1.0, (double)SPATIAL_ACC_MAX_CELLS_PER_DIM),
        (uint32_t)CLAMP(floor(norm_b / CELL_EXT), 1.0, (double)SPATIAL_ACC_MAX_CELLS_PER_DIM),
        (uint32_t)CLAMP(floor(norm_c / CELL_EXT), 1.0, (double)SPATIAL_ACC_MAX_CELLS_PER_DIM),
    };

#if DEBUG
    MD_LOG_DEBUG("cell_dim: %i %i %i", cell_dim[0], cell_dim[1], cell_dim[2]);
#endif

    store_frame(acc, A, I, origin, flags, cell_dim, inv_cell_ext, G00, G11, G22, H01, H02, H12);
    build_cells(acc, coords, in_flags);
}

void md_spatial_acc_init(md_spatial_acc_t* acc, const md_spatial_acc_desc_t* desc) {
    ASSERT(acc);
    ASSERT(desc);
    if (!desc->coords) {
        MD_LOG_ERROR("md_spatial_acc_init: the description has no coordinates");
        return;
    }
    spatial_acc_init_internal(acc, desc->coords, desc->cutoff, desc->unitcell, desc->flags);
}

// Generate forward neighbor offsets for a 3D grid cell
// - out: user-provided array of size at least (ncell*2+1)^3 - 1
// - ncell: number of cells in the neighborhood (1 → 1-cell, 2 → 2-cell, etc.)
// - Each element in out is int[4] representing offset {dx, dy, dz, 0}
// Returns the number of neighbors written
static inline size_t generate_forward_neighbors4(int out[][4], const int ncell[3]) {
    size_t count = 0;

    for (int dz = 0; dz <= ncell[2]; ++dz) {        // memory order: z-major
        for (int dy = -ncell[1]; dy <= ncell[1]; ++dy) {
            for (int dx = -ncell[0]; dx <= ncell[0]; ++dx) {
                // skip the origin
                if (dx == 0 && dy == 0 && dz == 0) continue;

                // Forward neighbor condition: i < j
                if (dz > 0 || (dz == 0 && dy > 0) || (dz == 0 && dy == 0 && dx > 0)) {
                    out[count][0] = dx;
                    out[count][1] = dy;
                    out[count][2] = dz;
					out[count][3] = 0;
                    count++;
                }
            }
        }
    }
    return count;
}

// Generate full neighbor offsets for a 3D grid cell
// - out: user-provided array of size at least (ncell*2+1)^3
// - ncell: number of cells in the neighborhood (1 → 1-cell, 2 → 2-cell, etc.)
// - Each element in out is int[4] representing offset {dx, dy, dz, 0}
// Returns the number of neighbors written
static inline size_t generate_neighbors4(int out[][4], const int ncell[3]) {
    size_t count = 0;
    for (int dz = -ncell[2]; dz <= ncell[2]; ++dz) {
        for (int dy = -ncell[1]; dy <= ncell[1]; ++dy) {
            for (int dx = -ncell[0]; dx <= ncell[0]; ++dx) {
                out[count][0] = dx;
                out[count][1] = dy;
                out[count][2] = dz;
                out[count][3] = 0;
                ++count;
            }
        }
    }
    return count;
}

static inline md_256 distance_squared_tri_256(md_256 dx, md_256 dy, md_256 dz, md_256 G00, md_256 G11, md_256 G22, md_256 H01, md_256 H02, md_256 H12) {
    md_256 dx2 = md_mm256_mul_ps(dx, dx);
    md_256 dy2 = md_mm256_mul_ps(dy, dy);
    md_256 dz2 = md_mm256_mul_ps(dz, dz);

    md_256 dxy = md_mm256_mul_ps(dx, dy);
    md_256 dxz = md_mm256_mul_ps(dx, dz);
    md_256 dyz = md_mm256_mul_ps(dy, dz);

    md_256 acc   = md_mm256_fmadd_ps(G00, dx2, md_mm256_fmadd_ps(G11, dy2, md_mm256_mul_ps(G22, dz2)));
    md_256 cross = md_mm256_fmadd_ps(H01, dxy, md_mm256_fmadd_ps(H02, dxz, md_mm256_mul_ps(H12, dyz)));
    return md_mm256_add_ps(acc, cross);
}

static inline md_256 distance_squared_ort_256(md_256 dx, md_256 dy, md_256 dz, md_256 G00, md_256 G11, md_256 G22) {
    md_256 dx2 = md_mm256_mul_ps(dx, dx);
    md_256 dy2 = md_mm256_mul_ps(dy, dy);
    md_256 dz2 = md_mm256_mul_ps(dz, dz);
    return md_mm256_fmadd_ps(G00, dx2, md_mm256_fmadd_ps(G11, dy2, md_mm256_mul_ps(G22, dz2)));
}

static inline int wrap_coord(int x, int N) {
    x += (x <  0) ? N : 0;
    x -= (x >= N) ? N : 0;
    return x;
}

static inline int isign(int a) {
    return (a > 0) - (a < 0);
}

// Number of cells the search has to reach along each axis, capped at ONE full period.
//
// The wrap applied to a neighbour index moves it by exactly one period, so it is only exact while
// the raw index (cell + offset) stays within [-N, 2N-1], that is while the offset magnitude stays
// within N. The cap is not an approximation: offsets spanning a full period in each direction
// already visit every cell of the grid at periodic shifts of -1, 0 and +1, which is the complete
// set of nearest image candidates. Reaching further only revisits those same cells in a more
// distant image, which cannot be nearer than one already considered.
//
// Without the cap a cutoff approaching the box size drives ncell past cell_dim. cell_dim is CLAMPed
// to at least 1, so it collapses precisely when the cutoff grows, and the single period wrap then
// hands back an index that is STILL out of range - which is then read straight into the element
// arrays.
static inline void neighbor_cell_extent(int out_ncell[3], double cutoff, const md_spatial_acc_t* acc) {
    for (int i = 0; i < 3; ++i) {
        const int dim = (int)acc->cell_dim[i];
        const double n = ceil(cutoff * (double)acc->inv_cell_ext[i] * (double)acc->cell_dim[i]);
        // Written so that a NaN lands on 0 rather than on an undefined float to int conversion.
        out_ncell[i] = !(n > 0.0) ? 0 : (n >= (double)dim ? dim : (int)n);
    }
}

static float calc_r2(double cutoff) {
    float r2 = (float)(cutoff * cutoff);
    return nextafterf(r2, r2 + 1.0f); // Round up to ensure we don't miss neighbors due to floating point precision
}

static inline vec4_t vec4_cart_to_fract(vec4_t in_c, const md_spatial_acc_t* acc) {
    const vec4_t t = vec4_set(acc->origin[0], acc->origin[1], acc->origin[2], 0);
    const vec4_t I[3] = {
        vec4_set(acc->I[0][0], acc->I[0][1], acc->I[0][2], 0),
        vec4_set(acc->I[1][0], acc->I[1][1], acc->I[1][2], 0),
        vec4_set(acc->I[2][0], acc->I[2][1], acc->I[2][2], 0),
    };
    return vec4_linear_combine_3(vec4_sub(in_c, t), I);
}

static inline void cart_to_fract(double out_s[3], const double in_c[3], const md_spatial_acc_t* acc) {
    // I is indexed as I[col][row]
    // fract = I * (cart - origin)
    const double x = in_c[0] - acc->origin[0];
    const double y = in_c[1] - acc->origin[1];
    const double z = in_c[2] - acc->origin[2];
    out_s[0] = acc->I[0][0] * x + acc->I[1][0] * y + acc->I[2][0] * z;
    out_s[1] = acc->I[0][1] * x + acc->I[1][1] * y + acc->I[2][1] * z;
    out_s[2] = acc->I[0][2] * x + acc->I[1][2] * y + acc->I[2][2] * z;
}

static inline void fract_to_cart(double out_x[3], const double in_s[3], const md_spatial_acc_t* acc) {
    // A is indexed as A[col][row]
    // cart = A * fract + origin
    out_x[0] = acc->A[0][0] * in_s[0] + acc->A[1][0] * in_s[1] + acc->A[2][0] * in_s[2] + acc->origin[0];
    out_x[1] = acc->A[0][1] * in_s[0] + acc->A[1][1] * in_s[1] + acc->A[2][1] * in_s[2] + acc->origin[1];
    out_x[2] = acc->A[0][2] * in_s[0] + acc->A[1][2] * in_s[1] + acc->A[2][2] * in_s[2] + acc->origin[2];
}

static inline void fract_to_cart_ort_256(
    md_256* out_cx, md_256* out_cy, md_256* out_cz,
    const md_256 in_sx, const md_256 in_sy, const md_256 in_sz,
    const md_256 A00, const md_256 A11, const md_256 A22,
    const md_256 O0,  const md_256 O1,  const md_256 O2)
{
    *out_cx = md_mm256_fmadd_ps(in_sx, A00, O0);
    *out_cy = md_mm256_fmadd_ps(in_sy, A11, O1);
    *out_cz = md_mm256_fmadd_ps(in_sz, A22, O2);
}

static inline void fract_to_cart_tri_256(
    md_256* out_cx, md_256* out_cy, md_256* out_cz,
    const md_256 in_sx, const md_256 in_sy, const md_256 in_sz,
    const md_256 A00, const md_256 A10, const md_256 A11, const md_256 A20, const md_256 A21, const md_256 A22,
    const md_256 O0,  const md_256 O1,  const md_256 O2)
{
    *out_cx = md_mm256_fmadd_ps(in_sx, A00, md_mm256_fmadd_ps(in_sy, A10, md_mm256_fmadd_ps(in_sz, A20, O0)));
    *out_cy = md_mm256_fmadd_ps(in_sy, A11, md_mm256_fmadd_ps(in_sz, A21, O1));
    *out_cz = md_mm256_fmadd_ps(in_sz, A22, O2);
}

static inline void batch_fract_to_cart_ort_256(float* x, float* y, float* z, size_t count, const md_spatial_acc_t* acc) {
	const md_256 A00 = md_mm256_set1_ps(acc->A[0][0]);
	const md_256 A11 = md_mm256_set1_ps(acc->A[1][1]);
	const md_256 A22 = md_mm256_set1_ps(acc->A[2][2]);

	const md_256 O0  = md_mm256_set1_ps(acc->origin[0]);
	const md_256 O1  = md_mm256_set1_ps(acc->origin[1]);
	const md_256 O2  = md_mm256_set1_ps(acc->origin[2]);

    for (size_t i = 0; i < count; i += 8) {
		md_256 sx = md_mm256_loadu_ps(x + i);
        md_256 sy = md_mm256_loadu_ps(y + i);
        md_256 sz = md_mm256_loadu_ps(z + i);

		md_256 cx, cy, cz;
        fract_to_cart_ort_256(&cx, &cy, &cz, sx, sy, sz, A00, A11, A22, O0, O1, O2);

		md_mm256_storeu_ps(x + i, cx);
		md_mm256_storeu_ps(y + i, cy);
		md_mm256_storeu_ps(z + i, cz);
    }
}

// The image of a query centre which the AABB query works in: folded into the cell along each periodic
// axis, left alone along the others. Both the bounds test and the returned coordinates are relative to it.
static inline void aabb_query_center(double out_frac[3], double out_cart[3], const double center[3], const md_spatial_acc_t* acc) {
    double s[3];
    cart_to_fract(s, center, acc);

    if (acc->flags & MD_UNITCELL_PBC_X) s[0] = fract(s[0]);
    if (acc->flags & MD_UNITCELL_PBC_Y) s[1] = fract(s[1]);
    if (acc->flags & MD_UNITCELL_PBC_Z) s[2] = fract(s[2]);

    if (out_frac) {
        out_frac[0] = s[0];
        out_frac[1] = s[1];
        out_frac[2] = s[2];
    }
    if (out_cart) {
        fract_to_cart(out_cart, s, acc);
    }
}

void md_spatial_acc_aabb_query_center(double out_center[3], const md_spatial_acc_t* acc, const double center[3]) {
    ASSERT(out_center);
    ASSERT(acc);
    ASSERT(center);
    aabb_query_center(NULL, out_center, center, acc);
}

static inline void cell_range_from_aabb_center_radius(
    int out_cmin[3],                 // inclusive
    int out_cmax[3],                 // exclusive (loop while ic < cmax)
    double out_fcen[3],              // fractional center (wrapped to [0,1) on periodic axes)
    double out_frad[3],              // fractional half-extents
    const double center[3],          // cartesian center
    const double radius[3],          // cartesian half-extents
    const md_spatial_acc_t* acc)
{
    const int pbc[3] = {
        (acc->flags & MD_UNITCELL_PBC_X) != 0,
        (acc->flags & MD_UNITCELL_PBC_Y) != 0,
        (acc->flags & MD_UNITCELL_PBC_Z) != 0,
    };
    
    // The same image md_spatial_acc_aabb_query_center reports - they have to agree exactly
    double sc[3];
    double cc[3];
    aabb_query_center(sc, cc, center, acc);

    // 2) Convert 8 corners to fractional, unwrap them near the center image, take component-wise min/max.
    double fmin[3] = { +DBL_MAX, +DBL_MAX, +DBL_MAX };
    double fmax[3] = { -DBL_MAX, -DBL_MAX, -DBL_MAX };

    for (int iz = 0; iz < 2; ++iz) {
        const double z = cc[2] + (iz ? +radius[2] : -radius[2]);
        for (int iy = 0; iy < 2; ++iy) {
            const double y = cc[1] + (iy ? +radius[1] : -radius[1]);
            for (int ix = 0; ix < 2; ++ix) {
                const double x = cc[0] + (ix ? +radius[0] : -radius[0]);
                const double c[3] = { x, y, z };

                double s[3];
                cart_to_fract(s, c, acc);

                for (int a = 0; a < 3; ++a) {
                    fmin[a] = MIN(fmin[a], s[a]);
                    fmax[a] = MAX(fmax[a], s[a]);
                }
            }
        }
    }

    out_fcen[0] = sc[0];
    out_fcen[1] = sc[1];
    out_fcen[2] = sc[2];

    out_frad[0] = 0.5 * (fmax[0] - fmin[0]);
    out_frad[1] = 0.5 * (fmax[1] - fmin[1]);
    out_frad[2] = 0.5 * (fmax[2] - fmin[2]);

    // 3) Convert fractional bounds to cell index bounds in the (fractional) grid.
    for (int a = 0; a < 3; ++a) {
        const int dim = (int)acc->cell_dim[a];

        int cmin = (int)floor(fmin[a] * (double)dim);
        int cmax = (int)ceil (fmax[a] * (double)dim);

        if (cmax <= cmin) cmax = cmin + 1;  // ensure at least one cell

        if (!pbc[a]) {
            cmin = CLAMP(cmin, 0, dim);
            cmax = CLAMP(cmax, 0, dim);
            if (cmax <= cmin) cmax = MIN(cmin + 1, dim);
        } else {
            // Cap the span at one full period either side of the centre cell. wrap_coord() moves an
            // index by exactly one period and the shift derived from it is only ever -1, 0 or +1, so
            // an index outside [-dim, 2*dim-1] comes back still out of range and is then read
            // straight into the element arrays. A half extent larger than the cell drives it there.
            //
            // Capping loses nothing: a full period in each direction already visits every cell of
            // the grid in each of its nearest images, and a more distant image cannot be nearer.
            const int ccen = (int)floor(out_fcen[a] * (double)dim);
            cmin = MAX(cmin, ccen - dim);
            cmax = MIN(cmax, ccen + dim + 1);
            if (cmax <= cmin) cmax = cmin + 1;
            ASSERT(cmin >= -dim && cmax <= 2 * dim);
        }

        out_cmin[a] = cmin;
        out_cmax[a] = cmax;
    }
}

// ### QUERY KERNELS ###
//
// Every pair query comes down to one operation: a single point, in fractional coordinates and moved to the periodic
// image the pair of cells needs, against a contiguous range of elements, eight at a time. The queries only differ in
// which ranges they pair up. Each kernel is written once and instantiated per metric (orthorhombic or triclinic)
// through a compile time constant argument, so each instantiation is specialized code with no branching on it at run
// time.

#define SPATIAL_ACC_BUFLEN 1024

// Results are staged here and handed to the callback in batches
typedef struct pair_out_t {
    md_spatial_acc_pair_callback_t callback;
    void*    user_param;
    size_t   count;
    uint32_t i[SPATIAL_ACC_BUFLEN];
    uint32_t j[SPATIAL_ACC_BUFLEN];
    float    d2[SPATIAL_ACC_BUFLEN];
} pair_out_t;

static inline void pair_out_flush(pair_out_t* out) {
    if (out->count) {
        out->callback(out->i, out->j, out->d2, out->count, out->user_param);
        out->count = 0;
    }
}

// Squared distance in the fractional frame and the cutoff, broadcast
typedef struct metric_t {
    md_256 G00, G11, G22;
    md_256 H01, H02, H12;
    md_256 r2;
} metric_t;

static inline metric_t metric_init(const md_spatial_acc_t* acc, double cutoff) {
    metric_t m;
    m.G00 = md_mm256_set1_ps(acc->G00);
    m.G11 = md_mm256_set1_ps(acc->G11);
    m.G22 = md_mm256_set1_ps(acc->G22);
    m.H01 = md_mm256_set1_ps(acc->H01);
    m.H02 = md_mm256_set1_ps(acc->H02);
    m.H12 = md_mm256_set1_ps(acc->H12);
    m.r2  = md_mm256_set1_ps(calc_r2(cutoff));
    return m;
}

// The element arrays, read once per query into locals. Read through acc instead, they would be reloaded after every
// store into the staging buffer, which as far as the compiler knows may alias anything in memory.
typedef struct elems_t {
    const float*    x;
    const float*    y;
    const float*    z;
    const uint32_t* idx;
} elems_t;

static inline elems_t elems_of(const md_spatial_acc_t* acc) {
    elems_t e = { acc->elem_x, acc->elem_y, acc->elem_z, acc->elem_idx };
    return e;
}

// Eight elements from ei/ex/ey/ez against the point (x, y, z): the pairs within the cutoff among the lanes set in
// lanes are staged. Takes the number of staged pairs and returns it updated, which keeps it in a register.
static FORCE_INLINE size_t point_vs_8(pair_out_t* out, size_t count, const metric_t* m, const float* ex, const float* ey, const float* ez, const uint32_t* ei,
                                      md_256 v_x, md_256 v_y, md_256 v_z, md_256i v_idx, int lanes, const bool tri) {
    const md_256 dx = md_mm256_sub_ps(v_x, md_mm256_loadu_ps(ex));
    const md_256 dy = md_mm256_sub_ps(v_y, md_mm256_loadu_ps(ey));
    const md_256 dz = md_mm256_sub_ps(v_z, md_mm256_loadu_ps(ez));
    md_256 d2 = tri ? distance_squared_tri_256(dx, dy, dz, m->G00, m->G11, m->G22, m->H01, m->H02, m->H12)
                    : distance_squared_ort_256(dx, dy, dz, m->G00, m->G11, m->G22);

    const int mask = md_mm256_movemask_ps(md_mm256_cmple_ps(d2, m->r2)) & lanes;
    if (mask) {
        if (count + 8 > SPATIAL_ACC_BUFLEN) {
            out->count = count;
            pair_out_flush(out);
            count = 0;
        }
        const md_256i perm = md_mm256_compression_mask_8x32(mask);
        md_mm256_storeu_epi32(out->i  + count, v_idx);
        md_mm256_storeu_epi32(out->j  + count, md_mm256_permutevar8x32_epi32(md_mm256_loadu_si256((const md_256i*)ei), perm));
        md_mm256_storeu_ps   (out->d2 + count, md_mm256_permutevar8x32_ps(d2, perm));
        count += popcnt32(mask);
    }
    return count;
}

// The point (x, y, z) with index idx against the elements [beg, beg + len): whole blocks of eight, then the rest
// with the lanes past the end masked off. The loads of the last block run up to 7 elements past the range, which the
// padding of the element arrays keeps in bounds.
static FORCE_INLINE size_t point_vs_range(pair_out_t* out, size_t count, const metric_t* m, elems_t e, md_256 v_x, md_256 v_y, md_256 v_z, md_256i v_idx, uint32_t beg, uint32_t len, const bool tri) {
    const float*    ex = e.x + beg;
    const float*    ey = e.y + beg;
    const float*    ez = e.z + beg;
    const uint32_t* ei = e.idx + beg;
    uint32_t j = 0;
    for (; j + 8 <= len; j += 8) {
        count = point_vs_8(out, count, m, ex + j, ey + j, ez + j, ei + j, v_x, v_y, v_z, v_idx, 0xFF, tri);
    }
    if (j < len) {
        count = point_vs_8(out, count, m, ex + j, ey + j, ez + j, ei + j, v_x, v_y, v_z, v_idx, (1 << (len - j)) - 1, tri);
    }
    return count;
}

// The elements of one cell against each other, every pair once
static FORCE_INLINE size_t cell_self(pair_out_t* out, size_t count, const metric_t* m, elems_t e, uint32_t off, uint32_t len, const bool tri) {
    for (uint32_t i = 0; i + 1 < len; ++i) {
        const uint32_t k = off + i;
        count = point_vs_range(out, count, m, e, md_mm256_set1_ps(e.x[k]), md_mm256_set1_ps(e.y[k]), md_mm256_set1_ps(e.z[k]),
                               md_mm256_set1_epi32((int)e.idx[k]), k + 1, len - i - 1, tri);
    }
    return count;
}

// The elements of cell i, moved by shift (in periods), against those of cell j
static FORCE_INLINE size_t cell_vs_cell(pair_out_t* out, size_t count, const metric_t* m, elems_t e, uint32_t off_i, uint32_t len_i, uint32_t off_j, uint32_t len_j, vec4_t shift, const bool tri) {
    const md_256 sx = md_mm256_set1_ps(shift.x);
    const md_256 sy = md_mm256_set1_ps(shift.y);
    const md_256 sz = md_mm256_set1_ps(shift.z);
    for (uint32_t i = 0; i < len_i; ++i) {
        const uint32_t k = off_i + i;
        count = point_vs_range(out, count, m, e,
                               md_mm256_add_ps(md_mm256_set1_ps(e.x[k]), sx),
                               md_mm256_add_ps(md_mm256_set1_ps(e.y[k]), sy),
                               md_mm256_add_ps(md_mm256_set1_ps(e.z[k]), sz),
                               md_mm256_set1_epi32((int)e.idx[k]), off_j, len_j, tri);
    }
    return count;
}

// Wrapping of neighbour cell coordinates into the grid
typedef struct grid_wrap_t {
    ivec4_t dim;
    ivec4_t dim_1;
    ivec4_t pbc;    // All bits set along the periodic axes
} grid_wrap_t;

static inline grid_wrap_t grid_wrap_init(const md_spatial_acc_t* acc) {
    grid_wrap_t g;
    g.dim   = ivec4_set((int)acc->cell_dim[0], (int)acc->cell_dim[1], (int)acc->cell_dim[2], 0);
    g.dim_1 = ivec4_sub(g.dim, ivec4_set(1, 1, 1, 0));
    g.pbc   = ivec4_set((acc->flags & MD_UNITCELL_PBC_X) ? -1 : 0, (acc->flags & MD_UNITCELL_PBC_Y) ? -1 : 0, (acc->flags & MD_UNITCELL_PBC_Z) ? -1 : 0, 0);
    return g;
}

// Cell coordinates n, at most one period outside the grid, wrapped into it. False when n lies beyond an edge which is
// not periodic. The shift is the image offset, in periods, to apply to the other side of the pair: +1 where n wrapped
// up from below, -1 where it wrapped down from above. It stays an integer vector: most neighbour cells of a sparse
// system turn out empty, and converting it is left to those which do not.
static FORCE_INLINE bool grid_wrap(const grid_wrap_t* g, ivec4_t n, uint32_t out_c[3], ivec4_t* out_shift) {
    const ivec4_t upper = ivec4_cmpgt(n, g->dim_1);
    const ivec4_t lower = ivec4_cmplt(n, ivec4_set1(0));
    if (ivec4_any(ivec4_andnot(ivec4_or(upper, lower), g->pbc))) return false;

    n = ivec4_add(n, ivec4_and(lower, g->dim));
    n = ivec4_sub(n, ivec4_and(upper, g->dim));
    int c[4];
    ivec4_store(c, n);
    out_c[0] = (uint32_t)c[0];
    out_c[1] = (uint32_t)c[1];
    out_c[2] = (uint32_t)c[2];
    *out_shift = ivec4_sub(ivec4_and(lower, ivec4_set1(1)), ivec4_and(upper, ivec4_set1(1)));
    return true;
}

// The elements of a cell, returned by value: through pointers to locals the compiler kept the length on the stack,
// and reloaded it in every iteration of the pair loop
typedef struct cell_range_t {
    uint32_t off;
    uint32_t len;
} cell_range_t;

// Recently used tier 1 nodes, for queries which look up cells one at a time near each other: the neighbours of an
// external point (which come in the order of the stream, mostly coherent), the cells of an AABB. Consecutive lookups
// fall in a few nodes, which then cost a hash and a compare instead of a walk from the top.
#define NODE_CACHE_SIZE 32

typedef struct node_cache_t {
    uint64_t key[NODE_CACHE_SIZE];      // Node coordinates, packed. UINT64_MAX when empty.
    uint64_t mask[NODE_CACHE_SIZE];
    uint32_t base[NODE_CACHE_SIZE];
} node_cache_t;

static inline void node_cache_init(node_cache_t* cache) {
    MEMSET(cache->key, 0xFF, sizeof(cache->key));
}

// The elements of the cell at c (within the grid), an empty range for an empty cell
static FORCE_INLINE cell_range_t cell_at(const md_spatial_acc_t* acc, node_cache_t* cache, const uint32_t c[3]) {
    const uint32_t nx = c[0] >> 2;
    const uint32_t ny = c[1] >> 2;
    const uint32_t nz = c[2] >> 2;
    const uint64_t key  = ((uint64_t)nz << 42) | ((uint64_t)ny << 21) | (uint64_t)nx;
    const uint32_t slot = ((nx * 73856093u) ^ (ny * 19349663u) ^ (nz * 83492791u)) & (NODE_CACHE_SIZE - 1);
    if (cache->key[slot] != key) {
        cache->key[slot] = key;
        node_lookup(acc, c, 1, &cache->mask[slot], &cache->base[slot]);
    }
    const uint64_t mask = cache->mask[slot];
    const uint32_t bit  = local_bit(c, 0);
    if (!((mask >> bit) & 1)) return (cell_range_t){ 0, 0 };
    const uint32_t k   = cache->base[slot] + (uint32_t)popcnt64(mask & ((1ULL << bit) - 1));
    const uint32_t off = acc->cell_off[k];
    return (cell_range_t){ off, acc->cell_off[k + 1] - off };
}

static bool neighbor_stencil_fits(const int ncell[3], const char* caller) {
    if (2 * ncell[0] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[1] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[2] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS) {
        MD_LOG_ERROR("%s: the cutoff is more than twice the cutoff the structure was built for", caller);
        return false;
    }
    return true;
}

// --- INTERNAL PAIRS ---

typedef struct cell_iter_t {
    const md_spatial_acc_t* acc;
    size_t   top;                                       // Next top node to visit
    uint32_t depth;                                     // Tier of the node whose children are visited, num_tiers + 1 between top nodes
    uint64_t rem[MD_SPATIAL_ACC_MAX_TIERS + 1];         // Children not yet visited, per tier
    uint32_t next[MD_SPATIAL_ACC_MAX_TIERS + 1];        // Index of the next child, per tier
    uint32_t org[MD_SPATIAL_ACC_MAX_TIERS + 1][3];      // Cell coordinates of the node's corner, per tier
    uint64_t mask[MD_SPATIAL_ACC_MAX_TIERS + 1];        // Mask and first child of the node, per tier
    uint32_t first[MD_SPATIAL_ACC_MAX_TIERS + 1];
} cell_iter_t;

static inline cell_iter_t cell_iter(const md_spatial_acc_t* acc) {
    cell_iter_t it;
    MEMSET(&it, 0, sizeof(it));
    it.acc = acc;
    it.depth = acc->num_tiers + 1;
    return it;
}

static inline bool cell_iter_next(cell_iter_t* it, uint32_t* out_idx, uint32_t out_c[3]) {
    const md_spatial_acc_t* acc = it->acc;
    const uint32_t L = acc->num_tiers;
    const size_t num_top = (size_t)acc->top_dim[0] * acc->top_dim[1] * acc->top_dim[2];
    for (;;) {
        const uint32_t d = it->depth;
        if (d > L) {
            while (it->top < num_top && acc->top_mask[it->top] == 0) ++it->top;
            if (it->top >= num_top) return false;
            const size_t t = it->top++;
            it->rem[L]  = it->mask[L]  = acc->top_mask[t];
            it->next[L] = it->first[L] = acc->top_base[t];
            const size_t txy = (size_t)acc->top_dim[0] * acc->top_dim[1];
            it->org[L][0] = (uint32_t)(t % acc->top_dim[0])       << (2 * L);
            it->org[L][1] = (uint32_t)((t % txy) / acc->top_dim[0]) << (2 * L);
            it->org[L][2] = (uint32_t)(t / txy)                   << (2 * L);
            it->depth = L;
            continue;
        }
        if (it->rem[d] == 0) {
            it->depth = d + 1;
            continue;
        }
        const uint32_t bit = (uint32_t)ctz64(it->rem[d]);
        it->rem[d] &= it->rem[d] - 1;
        const uint32_t idx = it->next[d]++;
        const uint32_t s = 2 * (d - 1);
        const uint32_t c[3] = {
            it->org[d][0] + ((bit & 3) << s),
            it->org[d][1] + (((bit >> 2) & 3) << s),
            it->org[d][2] + (((bit >> 4) & 3) << s),
        };
        if (d == 1) {
            *out_idx = idx;
            out_c[0] = c[0];
            out_c[1] = c[1];
            out_c[2] = c[2];
            return true;
        }
        it->rem[d - 1]  = it->mask[d - 1]  = acc->tier_mask[d - 1][idx];
        it->next[d - 1] = it->first[d - 1] = acc->tier_base[d - 1][idx];
        it->org[d - 1][0] = c[0];
        it->org[d - 1][1] = c[1];
        it->org[d - 1][2] = c[2];
        it->depth = d - 1;
    }
}

// Index of the occupied cell at cell coordinates c, or UINT32_MAX if it is empty. Starts from the lowest node on the
// iterator's current path which contains c rather than from the top: a neighbour of the current cell mostly shares
// its tier 1 or tier 2 node, which leaves one or two steps.
static inline uint32_t cell_lookup_near(const cell_iter_t* it, const uint32_t c[3]) {
    const md_spatial_acc_t* acc = it->acc;
    const uint32_t L = acc->num_tiers;
    uint32_t l = 1;
    for (; l <= L; ++l) {
        const uint32_t ext = 1u << (2 * l);
        // Unsigned: a coordinate below the corner wraps to a large value
        if (c[0] - it->org[l][0] < ext && c[1] - it->org[l][1] < ext && c[2] - it->org[l][2] < ext) break;
    }
    if (l > L) return cell_lookup(acc, c);
    uint64_t mask = it->mask[l];
    uint32_t base = it->first[l];
    for (uint32_t t = l - 1; ; --t) {
        const uint32_t bit = local_bit(c, t);
        if (!((mask >> bit) & 1)) return UINT32_MAX;
        const uint32_t idx = base + (uint32_t)popcnt64(mask & ((1ULL << bit) - 1));
        if (t == 0) return idx;
        mask = acc->tier_mask[t][idx];
        base = acc->tier_base[t][idx];
    }
}

// The node of tier T (>= 1) containing cell coordinates c: its mask and first child, a mask of 0 if there is none.
// Starts, as cell_lookup_near, from the lowest node on the iterator's current path which contains c.
static inline void node_lookup_near(const cell_iter_t* it, const uint32_t c[3], uint32_t T, uint64_t* out_mask, uint32_t* out_base) {
    const md_spatial_acc_t* acc = it->acc;
    const uint32_t L = acc->num_tiers;
    uint32_t l = T;
    for (; l <= L; ++l) {
        const uint32_t ext = 1u << (2 * l);
        if (c[0] - it->org[l][0] < ext && c[1] - it->org[l][1] < ext && c[2] - it->org[l][2] < ext) break;
    }
    if (l > L) {
        node_lookup(acc, c, T, out_mask, out_base);
        return;
    }
    uint64_t mask = it->mask[l];
    uint32_t base = it->first[l];
    for (uint32_t t = l; t > T; --t) {
        const uint32_t bit = local_bit(c, t - 1);
        if (!((mask >> bit) & 1)) {
            *out_mask = 0;
            *out_base = 0;
            return;
        }
        const uint32_t idx = base + (uint32_t)popcnt64(mask & ((1ULL << bit) - 1));
        mask = acc->tier_mask[t - 1][idx];
        base = acc->tier_base[t - 1][idx];
    }
    *out_mask = mask;
    *out_base = base;
}


static FORCE_INLINE void internal_pairs(const md_spatial_acc_t* acc, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param, const bool tri) {
    int ncell[3];
    neighbor_cell_extent(ncell, cutoff, acc);
    if (!neighbor_stencil_fits(ncell, "md_spatial_acc_for_each_internal_pair_within_cutoff")) return;

    int nbr[SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS][4];
    const size_t num_nbr = generate_forward_neighbors4(nbr, ncell);

    const metric_t    m = metric_init(acc, cutoff);
    const grid_wrap_t g = grid_wrap_init(acc);
    const elems_t     e = elems_of(acc);
    pair_out_t out;
    out.callback   = callback;
    out.user_param = user_param;
    out.count      = 0;
    size_t count   = 0;

    // The tier 1 nodes around the current one (3x3x3), filled as the neighbours ask for them. A neighbour within the
    // grid is at most two cells away, so its tier 1 node (4 cells across) is one of these.
    uint64_t halo_mask[27];
    uint32_t halo_base[27];
    uint32_t halo_valid = 0;
    uint32_t cur_node[3] = { UINT32_MAX, UINT32_MAX, UINT32_MAX };

    cell_iter_t it = cell_iter(acc);
    uint32_t ci;
    uint32_t cc[3];
    while (cell_iter_next(&it, &ci, cc)) {
        if ((cc[0] >> 2) != cur_node[0] || (cc[1] >> 2) != cur_node[1] || (cc[2] >> 2) != cur_node[2]) {
            cur_node[0] = cc[0] >> 2;
            cur_node[1] = cc[1] >> 2;
            cur_node[2] = cc[2] >> 2;
            halo_valid = 0;
        }
        const uint32_t off_i = acc->cell_off[ci];
        const uint32_t len_i = acc->cell_off[ci + 1] - off_i;

        count = cell_self(&out, count, &m, e, off_i, len_i, tri);

        const ivec4_t c_v = ivec4_set((int)cc[0], (int)cc[1], (int)cc[2], 0);
        for (size_t n = 0; n < num_nbr; ++n) {
            const ivec4_t n_v = ivec4_add(c_v, ivec4_load(nbr[n]));
            uint32_t nc[3];
            ivec4_t shift;
            if (!grid_wrap(&g, n_v, nc, &shift)) continue;

            uint32_t cj;
            if (!ivec4_any(shift)) {
                // Within the grid: through the cached tier 1 nodes
                const uint32_t slot = ((nc[0] >> 2) - cur_node[0] + 1) + 3 * ((nc[1] >> 2) - cur_node[1] + 1) + 9 * ((nc[2] >> 2) - cur_node[2] + 1);
                ASSERT(slot < 27);
                if (!((halo_valid >> slot) & 1)) {
                    halo_valid |= 1u << slot;
                    node_lookup_near(&it, nc, 1, &halo_mask[slot], &halo_base[slot]);
                }
                const uint32_t bit = local_bit(nc, 0);
                if (!((halo_mask[slot] >> bit) & 1)) continue;
                cj = halo_base[slot] + (uint32_t)popcnt64(halo_mask[slot] & ((1ULL << bit) - 1));
            } else {
                cj = cell_lookup_near(&it, nc);
                if (cj == UINT32_MAX) continue;
            }
            const uint32_t off_j = acc->cell_off[cj];
            const uint32_t len_j = acc->cell_off[cj + 1] - off_j;
            count = cell_vs_cell(&out, count, &m, e, off_i, len_i, off_j, len_j, vec4_from_ivec4(shift), tri);
        }
    }
    out.count = count;
    pair_out_flush(&out);
}

static void internal_pairs_ortho(const md_spatial_acc_t* acc, double cutoff, md_spatial_acc_pair_callback_t cb, void* user) { internal_pairs(acc, cutoff, cb, user, false); }
static void internal_pairs_tricl(const md_spatial_acc_t* acc, double cutoff, md_spatial_acc_pair_callback_t cb, void* user) { internal_pairs(acc, cutoff, cb, user, true);  }

// --- EXTERNAL PAIRS ---

static FORCE_INLINE void external_pairs(const md_spatial_acc_t* acc, const md_coord_stream_t* ext, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param, md_spatial_acc_flags_t flags, const bool tri) {
    int ncell[3];
    neighbor_cell_extent(ncell, cutoff, acc);
    if (!neighbor_stencil_fits(ncell, "md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff")) return;

    if ((flags & MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX) && !ext->idx) {
        MD_LOG_ERROR("md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff: MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX is set but the external stream has no idx");
        return;
    }

    int nbr[SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS][4];
    const size_t num_nbr = generate_neighbors4(nbr, ncell);

    const metric_t    m = metric_init(acc, cutoff);
    const grid_wrap_t g = grid_wrap_init(acc);
    const elems_t     e = elems_of(acc);
    pair_out_t out;
    out.callback   = callback;
    out.user_param = user_param;
    out.count      = 0;
    size_t count   = 0;

    node_cache_t cache;
    node_cache_init(&cache);

    vec4_t fract_mask;
    MEMCPY(&fract_mask, &g.pbc, sizeof(fract_mask));
    const vec4_t fcell_dim = vec4_set((float)acc->cell_dim[0], (float)acc->cell_dim[1], (float)acc->cell_dim[2], 0);

    for (size_t pt = 0; pt < ext->count; ++pt) {
        const uint32_t idx = (flags & MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX) ? (uint32_t)md_coord_stream_load_idx(ext, pt) : (uint32_t)pt;
        // Into the cell along the periodic axes. Along the others a point outside the grid has no cells to visit
        // beyond the edge, which the wrapping rejects.
        vec4_t f = vec4_cart_to_fract(md_coord_stream_load_vec4(ext, pt), acc);
        f = vec4_blend(f, vec4_fract(f), fract_mask);
        const ivec4_t c_v = ivec4_from_vec4(vec4_floor(vec4_mul(f, fcell_dim)));
        const md_256i v_idx = md_mm256_set1_epi32((int)idx);

        for (size_t n = 0; n < num_nbr; ++n) {
            uint32_t nc[3];
            ivec4_t shift;
            if (!grid_wrap(&g, ivec4_add(c_v, ivec4_load(nbr[n])), nc, &shift)) continue;
            const cell_range_t r = cell_at(acc, &cache, nc);
            if (!r.len) continue;
            const vec4_t fs = vec4_add(f, vec4_from_ivec4(shift));
            count = point_vs_range(&out, count, &m, e, md_mm256_set1_ps(fs.x), md_mm256_set1_ps(fs.y), md_mm256_set1_ps(fs.z), v_idx, r.off, r.len, tri);
        }
    }
    out.count = count;
    pair_out_flush(&out);
}

static void external_pairs_ortho(const md_spatial_acc_t* acc, const md_coord_stream_t* ext, double cutoff, md_spatial_acc_pair_callback_t cb, void* user, md_spatial_acc_flags_t flags) { external_pairs(acc, ext, cutoff, cb, user, flags, false); }
static void external_pairs_tricl(const md_spatial_acc_t* acc, const md_coord_stream_t* ext, double cutoff, md_spatial_acc_pair_callback_t cb, void* user, md_spatial_acc_flags_t flags) { external_pairs(acc, ext, cutoff, cb, user, flags, true);  }

// --- POINTS IN AABB ---

typedef struct point_out_t {
    md_spatial_acc_point_callback_t callback;
    void*    user_param;
    size_t   count;
    uint32_t i[SPATIAL_ACC_BUFLEN];
    float    x[SPATIAL_ACC_BUFLEN];
    float    y[SPATIAL_ACC_BUFLEN];
    float    z[SPATIAL_ACC_BUFLEN];
} point_out_t;

// The orthorhombic kernel stages fractional coordinates and converts them a batch at a time, the triclinic one has
// cartesian coordinates already (it tests the bounds in them)
static FORCE_INLINE void point_out_flush(point_out_t* out, const md_spatial_acc_t* acc, const bool tri) {
    if (out->count) {
        if (!tri) batch_fract_to_cart_ort_256(out->x, out->y, out->z, out->count, acc);
        out->callback(out->i, out->x, out->y, out->z, out->count, out->user_param);
        out->count = 0;
    }
}

static FORCE_INLINE void points_in_aabb(const md_spatial_acc_t* acc, const double aabb_cen[3], const double aabb_rad[3], md_spatial_acc_point_callback_t callback, void* user_param, const bool tri) {
    point_out_t out;
    out.callback   = callback;
    out.user_param = user_param;
    out.count      = 0;

    node_cache_t cache;
    node_cache_init(&cache);

    const uint32_t* cdim = acc->cell_dim;
    const int pbc[3] = {
        (acc->flags & MD_UNITCELL_PBC_X) != 0,
        (acc->flags & MD_UNITCELL_PBC_Y) != 0,
        (acc->flags & MD_UNITCELL_PBC_Z) != 0,
    };

    int cmin[3], cmax[3];
    double fcen[3], frad[3];
    cell_range_from_aabb_center_radius(cmin, cmax, fcen, frad, aabb_cen, aabb_rad, acc);

    // Orthorhombic: the box is a box in the fractional frame as well, tested there. Triclinic: it is not, so the
    // elements are taken to cartesian coordinates and tested against the box itself.
    md_256 v_min[3], v_max[3];
    if (!tri) {
        for (int a = 0; a < 3; ++a) {
            // Along a periodic axis half a period already reaches every image; along one which is not, the frame is
            // only the extent of the points, and the box may well reach past it
            const double r = pbc[a] ? MIN(frad[a], 0.5) : frad[a];
            v_min[a] = md_mm256_set1_ps((float)(fcen[a] - r));
            v_max[a] = md_mm256_set1_ps((float)(fcen[a] + r));
        }
    } else {
        double cc[3];
        fract_to_cart(cc, fcen, acc);
        for (int a = 0; a < 3; ++a) {
            v_min[a] = md_mm256_set1_ps((float)(cc[a] - aabb_rad[a]));
            v_max[a] = md_mm256_set1_ps((float)(cc[a] + aabb_rad[a]));
        }
    }
    const md_256 A00 = md_mm256_set1_ps(acc->A[0][0]);
    const md_256 A10 = md_mm256_set1_ps(acc->A[1][0]);
    const md_256 A11 = md_mm256_set1_ps(acc->A[1][1]);
    const md_256 A20 = md_mm256_set1_ps(acc->A[2][0]);
    const md_256 A21 = md_mm256_set1_ps(acc->A[2][1]);
    const md_256 A22 = md_mm256_set1_ps(acc->A[2][2]);
    const md_256 O0  = md_mm256_set1_ps(acc->origin[0]);
    const md_256 O1  = md_mm256_set1_ps(acc->origin[1]);
    const md_256 O2  = md_mm256_set1_ps(acc->origin[2]);

    const md_256i add8 = md_mm256_set1_epi32(8);

    for (int icz = cmin[2]; icz < cmax[2]; ++icz) {
        const int cz = pbc[2] ? wrap_coord(icz, (int)cdim[2]) : icz;
        const md_256 shift_z = md_mm256_set1_ps((float)isign(icz - cz));

        for (int icy = cmin[1]; icy < cmax[1]; ++icy) {
            const int cy = pbc[1] ? wrap_coord(icy, (int)cdim[1]) : icy;
            const md_256 shift_y = md_mm256_set1_ps((float)isign(icy - cy));

            for (int icx = cmin[0]; icx < cmax[0]; ++icx) {
                const int cx = pbc[0] ? wrap_coord(icx, (int)cdim[0]) : icx;
                const md_256 shift_x = md_mm256_set1_ps((float)isign(icx - cx));

                const uint32_t c[3] = { (uint32_t)cx, (uint32_t)cy, (uint32_t)cz };
                const cell_range_t r = cell_at(acc, &cache, c);
                if (!r.len) continue;
                const uint32_t off = r.off;
                const uint32_t len = r.len;

                const float*    ex = acc->elem_x + off;
                const float*    ey = acc->elem_y + off;
                const float*    ez = acc->elem_z + off;
                const uint32_t* ei = acc->elem_idx + off;
                const md_256i v_len = md_mm256_set1_epi32((int)len);
                md_256i v_j = md_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
                size_t count = out.count;

                for (uint32_t j = 0; j < len; j += 8) {
                    const md_256 j_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len));

                    // Fractional coordinates in the image of the visited cell
                    md_256 px = md_mm256_add_ps(md_mm256_loadu_ps(ex + j), shift_x);
                    md_256 py = md_mm256_add_ps(md_mm256_loadu_ps(ey + j), shift_y);
                    md_256 pz = md_mm256_add_ps(md_mm256_loadu_ps(ez + j), shift_z);
                    if (tri) {
                        md_256 qx, qy, qz;
                        fract_to_cart_tri_256(&qx, &qy, &qz, px, py, pz, A00, A10, A11, A20, A21, A22, O0, O1, O2);
                        px = qx;
                        py = qy;
                        pz = qz;
                    }

                    const md_256 in_min = md_mm256_and_ps(md_mm256_cmpge_ps(px, v_min[0]), md_mm256_and_ps(md_mm256_cmpge_ps(py, v_min[1]), md_mm256_cmpge_ps(pz, v_min[2])));
                    const md_256 in_max = md_mm256_and_ps(md_mm256_cmple_ps(px, v_max[0]), md_mm256_and_ps(md_mm256_cmple_ps(py, v_max[1]), md_mm256_cmple_ps(pz, v_max[2])));
                    const int mask = md_mm256_movemask_ps(md_mm256_and_ps(md_mm256_and_ps(in_min, in_max), j_mask));
                    if (mask) {
                        if (count + 8 > SPATIAL_ACC_BUFLEN) {
                            out.count = count;
                            point_out_flush(&out, acc, tri);
                            count = 0;
                        }

                        const md_256i perm = md_mm256_compression_mask_8x32(mask);
                        md_mm256_storeu_epi32(out.i + count, md_mm256_permutevar8x32_epi32(md_mm256_loadu_si256((const md_256i*)(ei + j)), perm));
                        md_mm256_storeu_ps   (out.x + count, md_mm256_permutevar8x32_ps(px, perm));
                        md_mm256_storeu_ps   (out.y + count, md_mm256_permutevar8x32_ps(py, perm));
                        md_mm256_storeu_ps   (out.z + count, md_mm256_permutevar8x32_ps(pz, perm));
                        count += popcnt32(mask);
                    }
                    v_j = md_mm256_add_epi32(v_j, add8);
                }
                out.count = count;
            }
        }
    }
    point_out_flush(&out, acc, tri);
}

static void points_in_aabb_ortho(const md_spatial_acc_t* acc, const double cen[3], const double rad[3], md_spatial_acc_point_callback_t cb, void* user) { points_in_aabb(acc, cen, rad, cb, user, false); }
static void points_in_aabb_tricl(const md_spatial_acc_t* acc, const double cen[3], const double rad[3], md_spatial_acc_point_callback_t cb, void* user) { points_in_aabb(acc, cen, rad, cb, user, true);  }

// ### PUBLIC QUERIES ###

void md_spatial_acc_for_each_internal_pair_within_cutoff(const md_spatial_acc_t* acc, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param) {
    ASSERT(acc);
    ASSERT(callback);
    if (acc->num_elems == 0) return;
    const bool tricl = (acc->flags & MD_UNITCELL_TRICLINIC) != 0;
    tricl ? internal_pairs_tricl(acc, cutoff, callback, user_param) : internal_pairs_ortho(acc, cutoff, callback, user_param);
}

void md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(const md_spatial_acc_t* acc, const md_coord_stream_t* ext_stream, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param, md_spatial_acc_flags_t flags) {
    ASSERT(acc);
    ASSERT(ext_stream);
    ASSERT(callback);
    if (acc->num_elems == 0) return;
    const bool tricl = (acc->flags & MD_UNITCELL_TRICLINIC) != 0;
    tricl ? external_pairs_tricl(acc, ext_stream, cutoff, callback, user_param, flags) : external_pairs_ortho(acc, ext_stream, cutoff, callback, user_param, flags);
}

void md_spatial_acc_for_each_point_in_aabb(const md_spatial_acc_t* acc, const double aabb_cen[3], const double aabb_rad[3], md_spatial_acc_point_callback_t callback, void* user_param) {
    ASSERT(acc);
    ASSERT(callback);
    for (int i = 0; i < 3; ++i) {
        if (aabb_rad[i] < 0) {
            MD_LOG_ERROR("md_spatial_acc_for_each_point_in_aabb: negative radius not allowed");
            return;
        }
    }
    if (acc->num_elems == 0) return;
    const bool tricl = (acc->flags & MD_UNITCELL_TRICLINIC) != 0;
    tricl ? points_in_aabb_tricl(acc, aabb_cen, aabb_rad, callback, user_param) : points_in_aabb_ortho(acc, aabb_cen, aabb_rad, callback, user_param);
}

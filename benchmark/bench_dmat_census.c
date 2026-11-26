// md_dmat_census_bench
//
// What a density matrix kept in memory at all times, in 16 x 16 tiles of the AO pairs that matter, would hold, and
// what leaving the other pairs out does: a census of the patterns each consumer needs, on the CPU (no GPU), on
// VeloxChem files.
//
//   load      the times to read the dense matrix and the SCF density's orbitals from the file
//   patterns  for each pattern its shell pairs, AO pairs and tile pairs (the tiles in the AO order of md_gto_gpu_dmat:
//             radius class at 1e-6, then atom, then shell), their bytes in float and in double (with the map),
//             and how much of viamd's GPU density pattern it holds:
//               geometric   shells whose centres are no further apart than their radii at a cutoff plus a margin
//                           (md_gto_gpu_dmat): the density's (1e-6, margin 0 and 1 bohr, the latter viamd's) and
//                           ELF/LOL's (1e-7, margin 1)
//               topology    what md_topo_compute_extremum_graph_gto reads (md_topo_gto_reach): its local lists
//                           (local radius + region margin), its polish (far radius), both, and its GPU sweep's tiles
//                           as made today (far radius + region margin)
//               D-weighted  a shell pair kept where max |D| over the pair times a bound on sup |phi_a phi_b| (radial
//                           envelopes r^l sum |c| exp(-a r^2), never below the functions) exceeds delta; with E, a
//                           bound on what all left-out pairs together add to rho anywhere, |rho - rho'| <= E (rho
//                           only: the derivatives' bounds are not taken)
//   ESP       the charge distribution (md_gto_int_charges_init, viamd's threshold 1e-6 and far field) from D with
//             only a pattern's shell pairs (the tiles hold more: this is the worst case), against the whole D: the
//             gaussians, the build's time, the net charge's change, the potential's largest and RMS change at
//             points near the surface (--probe bohr from the nearest atom), and for scale the far field's own
//             error at the first --exact points
//
//   md_dmat_census_bench [options] [file ...]
//
//   --deltas <list>  D-weighted thresholds, comma separated (default 1e-6,1e-8,1e-10,1e-12,1e-14,1e-16)
//   --points <n>     ESP points (default 256)
//   --exact <n>      of them evaluated without the far field too (default 16; 0 skips)
//   --probe <bohr>   the points' distance from the nearest atom (default 4)
//   --no-esp, --no-topo
//   --threads <n>    (default: the machine's)
//   --data <dir>     test_data directory (default: the source tree's)
//
// A file is a path relative to <data>/vlx or an existing file. Without any: protein.h5 if present, else mol.h5.

#include <md_vlx.h>
#include <md_system.h>
#include <md_gto.h>
#include <md_gto_int.h>
#include <md_qm.h>
#include <md_topo.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_unit.h>

#include "bench_rev.h"      // MD_BENCH_SOURCE_REV, generated at build time by bench_version.cmake
#include "../unittest/test_executor.h"
#include "../unittest/qm_test_util.h"   // qm_test_dmat_from_matrix: a whole matrix as tiles, which only a measurement needs

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Bumped with every change to this program's measurement or output
#define MD_DMAT_CENSUS_BENCH_REV "dmat-census r1"

#define ANGSTROM_TO_BOHR 1.8897261246257702
#define TILE 16u
#define MAX_DELTAS 12

#define ENV_H 0.25      // bohr: the envelopes' grid
#define ENV_K 513       // 0 .. 128 bohr

static double ms_since(md_tick_t t0) { return md_tick_to_milliseconds(md_tick_now() - t0); }

// ---------------------------------------------------------------------------------------------------- the input

typedef struct src_t {
    md_allocator_i* arena;
    md_gto_basis_t  basis;
    float*          xyz;        // bohr
    double*         Z;          // nuclear charges, e
    size_t          num_atoms, num_ao;
    double*         D;          // total density [num_ao][num_ao], heap
    size_t          num_occ;
    double          ms_read_d, ms_read_orbitals, ms_open;
} src_t;

static void src_free(src_t* s) {
    if (s->D) md_free(md_get_heap_allocator(), s->D, sizeof(double) * s->num_ao * s->num_ao);
    if (s->arena) md_arena_allocator_destroy(s->arena);
    memset(s, 0, sizeof(*s));
}

static bool src_load(src_t* s, const char* path) {
    memset(s, 0, sizeof(*s));
    s->arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(64));
    md_system_t sys = { .alloc = s->arena };
    md_system_state_t state = { .alloc = s->arena };
    md_tick_t t0 = md_tick_now();
    if (!md_vlx_system_init_from_file(&sys, &state, md_str_from_cstr(path), NULL)) return false;
    s->ms_open = ms_since(t0);
    const md_attribute_t* c = md_attributes_find(&sys.attributes, MD_STR_LIT("qm/atom/coordinate"));
    const md_attribute_t* z = md_attributes_find(&sys.attributes, MD_STR_LIT("qm/atom/nuclear_charge"));
    const md_attribute_t* d = md_attributes_find(&sys.attributes, MD_STR_LIT("orbital/total/density"));
    double scale = 1.0;
    if (!d) {
        d = md_attributes_find(&sys.attributes, MD_STR_LIT("orbital/alpha/density"));
        scale = 2.0;    // a closed shell's alpha density: half the total
    }
    if (!c || !z || !d || !md_gto_basis_extract_attributes(&s->basis, &sys.attributes, s->arena)) {
        fprintf(stderr, "%s: no basis, geometry, nuclear charges or density\n", path);
        return false;
    }
    s->num_atoms = md_attribute_value_count(&c->format);
    double* xyz = md_arena_allocator_push(s->arena, sizeof(double) * 3 * s->num_atoms);
    md_attribute_extract_f64(xyz, 3 * s->num_atoms, &sys.attributes, c, md_attribute_slice_all(), md_unit_none());
    s->xyz = md_arena_allocator_push(s->arena, sizeof(float) * 3 * s->num_atoms);
    for (size_t i = 0; i < 3 * s->num_atoms; ++i) s->xyz[i] = (float)(xyz[i] * ANGSTROM_TO_BOHR);
    if (md_attribute_value_count(&z->format) < s->num_atoms) return false;
    s->Z = md_arena_allocator_push(s->arena, sizeof(double) * s->num_atoms);
    md_attribute_extract_f64(s->Z, s->num_atoms, &sys.attributes, z, md_attribute_slice_all(), md_unit_elementary_charge());
    s->num_ao = md_gto_basis_num_ao(&s->basis);
    const size_t nd = s->num_ao * s->num_ao;
    s->D = md_alloc(md_get_heap_allocator(), sizeof(double) * nd);
    t0 = md_tick_now();
    const size_t got = md_attribute_extract_f64(s->D, nd, &sys.attributes, d, md_attribute_slice_all(), md_unit_none());
    s->ms_read_d = ms_since(t0);
    if (got != nd) {
        fprintf(stderr, "%s: the density is not %zu x %zu\n", path, s->num_ao, s->num_ao);
        return false;
    }
    if (scale != 1.0) for (size_t i = 0; i < nd; ++i) s->D[i] *= scale;
    // the SCF density's orbitals, as viamd builds its GPU tiles from them: only the read is timed
    double* rows = NULL;
    double* w = NULL;
    size_t nao = 0;
    t0 = md_tick_now();
    s->num_occ = md_qm_density_orbitals(&rows, &w, &nao, &sys, MD_STR_LIT("orbital/total/density"), md_get_heap_allocator());
    s->ms_read_orbitals = ms_since(t0);
    if (rows) md_free(md_get_heap_allocator(), rows, sizeof(double) * s->num_occ * nao);
    if (w) md_free(md_get_heap_allocator(), w, sizeof(double) * s->num_occ);
    return true;
}

// ---------------------------------------------------------------------------------------------------- shells

typedef struct shells_t {
    uint32_t  ns, nao, nt;
    uint32_t* atom;
    uint32_t* l;
    uint32_t* nc;
    uint32_t* ao;       // first AO (the basis' order)
    uint32_t* tao0;     // first AO in md_gto_gpu_dmat's tiled order
    double*   A;        // 3 per shell, bohr
} shells_t;

// A shell's screening radius at a cutoff, exactly as md_gto_gpu_basis and md_gto_gpu_dmat make it
// (gto_expand_basis_gpu_shells: the largest over components and primitives, the coefficient rounded to float)
static float shell_radius(const md_gto_basis_t* b, uint32_t s, double cutoff) {
    const md_gto_shell_t* sh = &b->shells[s];
    const int l = (int)sh->l;
    double max_r = 0.0;
    for (int i = l; i >= 0; --i) for (int j = l - i; j >= 0; --j) {
        const int k = l - i - j;
        const double nrm = md_gto_cart_norm_factor(i, j, k);
        for (uint32_t p = 0; p < sh->num_primitives; ++p) {
            const float alpha = b->alpha[sh->primitive_offset + p];
            const float coeff = (float)(b->coeff[sh->primitive_offset + p] * nrm);
            const double r = md_gto_compute_radius_of_influence(i, j, k, (double)coeff, (double)alpha, cutoff);
            if (r > max_r) max_r = r;
        }
    }
    return (float)max_r;
}

static int cmp_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return x < y ? -1 : x > y;
}

// The tiled order: md_gto_gpu_dmat's at order_cutoff (radius class of 2 bohr, then atom, then shell)
static void shells_make(shells_t* S, const src_t* src, double order_cutoff) {
    const md_gto_basis_t* b = &src->basis;
    const uint32_t ns = b->num_shells;
    S->ns = ns;
    S->nao = (uint32_t)src->num_ao;
    S->nt = (S->nao + TILE - 1) / TILE;
    S->atom = calloc(ns, sizeof(uint32_t));
    S->l = calloc(ns, sizeof(uint32_t));
    S->nc = calloc(ns, sizeof(uint32_t));
    S->ao = calloc(ns, sizeof(uint32_t));
    S->tao0 = calloc(ns, sizeof(uint32_t));
    S->A = calloc(3 * (size_t)ns, sizeof(double));
    uint64_t* key = calloc(ns, sizeof(uint64_t));
    uint32_t acc = 0;
    for (uint32_t s = 0; s < ns; ++s) {
        S->atom[s] = b->shells[s].atom_idx;
        S->l[s] = b->shells[s].l;
        S->nc[s] = md_gto_num_cart_ao(b->shells[s].l);
        S->ao[s] = acc;
        acc += S->nc[s];
        for (int k = 0; k < 3; ++k) S->A[3 * s + k] = src->xyz[3 * S->atom[s] + k];
        const float r = shell_radius(b, s, order_cutoff);
        key[s] = ((uint64_t)(r * 0.5f) << 40) | ((uint64_t)S->atom[s] << 20) | s;
    }
    qsort(key, ns, sizeof(uint64_t), cmp_u64);
    acc = 0;
    for (uint32_t i = 0; i < ns; ++i) {
        const uint32_t s = (uint32_t)(key[i] & 0xFFFFFu);
        S->tao0[s] = acc;
        acc += S->nc[s];
    }
    free(key);
}

static void shells_free(shells_t* S) {
    free(S->atom); free(S->l); free(S->nc); free(S->ao); free(S->tao0); free(S->A);
    memset(S, 0, sizeof(*S));
}

static double shell_dist2(const shells_t* S, uint32_t s, uint32_t t) {
    const double dx = S->A[3 * s] - S->A[3 * t], dy = S->A[3 * s + 1] - S->A[3 * t + 1], dz = S->A[3 * s + 2] - S->A[3 * t + 2];
    return dx * dx + dy * dy + dz * dz;
}

// ---------------------------------------------------------------------------------------------------- patterns
// A set of shell pairs s <= t: one bit each, rows padded to whole words so that rows can be written in parallel.

typedef struct pattern_t {
    char      name[96];
    uint32_t  ns, words;    // words per row
    uint64_t* bits;
} pattern_t;

static void pattern_init(pattern_t* p, const char* name, uint32_t ns) {
    memset(p, 0, sizeof(*p));
    snprintf(p->name, sizeof(p->name), "%s", name);
    p->ns = ns;
    p->words = (ns + 63) / 64;
    p->bits = calloc((size_t)ns * p->words, sizeof(uint64_t));
}
static void pattern_free(pattern_t* p) { free(p->bits); memset(p, 0, sizeof(*p)); }
static inline void pattern_set(pattern_t* p, uint32_t s, uint32_t t) {
    if (s > t) { const uint32_t u = s; s = t; t = u; }
    p->bits[(size_t)s * p->words + t / 64] |= 1ull << (t % 64);
}
static inline bool pattern_has(const pattern_t* p, uint32_t s, uint32_t t) {
    if (s > t) { const uint32_t u = s; s = t; t = u; }
    return (p->bits[(size_t)s * p->words + t / 64] >> (t % 64)) & 1u;
}

// Shells whose centres are no further apart than ra + rb + margin (and both radii > 0, as md_gto_gpu_dmat_create)
static void pattern_geometric(pattern_t* p, const shells_t* S, const double* ra, const double* rb_or_null, double margin) {
    const double* rb = rb_or_null ? rb_or_null : ra;
    for (uint32_t s = 0; s < S->ns; ++s) {
        if (!(ra[s] > 0.0)) continue;
        for (uint32_t t = s; t < S->ns; ++t) {
            if (!(rb[t] > 0.0)) continue;
            const double lim = ra[s] + rb[t] + margin;
            if (shell_dist2(S, s, t) <= lim * lim) pattern_set(p, s, t);
        }
    }
}

typedef struct census_t {
    uint64_t shell_pairs, ao_pairs, tiles;
    double   ref_held;      // of the reference pattern's AO pairs, the fraction held
    double   E;             // the bound on what the pairs left out add to rho anywhere (wsum summed), < 0 without
} census_t;

static census_t pattern_census(const pattern_t* p, const shells_t* S, const pattern_t* ref, const uint64_t* row, const double* wsum) {
    census_t c = {0};
    c.E = wsum ? 0.0 : -1.0;
    uint8_t* on = calloc((size_t)S->nt * S->nt, 1);
    uint64_t ref_pairs = 0, ref_in = 0;
    for (uint32_t s = 0; s < S->ns; ++s) {
        for (uint32_t t = s; t < S->ns; ++t) {
            const uint64_t pairs = s == t ? (uint64_t)S->nc[s] * (S->nc[s] + 1) / 2 : (uint64_t)S->nc[s] * S->nc[t];
            const bool has = pattern_has(p, s, t);
            if (ref && pattern_has(ref, s, t)) { ref_pairs += pairs; if (has) ref_in += pairs; }
            if (!has) {
                if (wsum && t > s) c.E += wsum[row[s] + (t - s - 1)];
                continue;
            }
            c.shell_pairs++;
            c.ao_pairs += pairs;
            for (uint32_t i = S->tao0[s] / TILE; i <= (S->tao0[s] + S->nc[s] - 1) / TILE; ++i) {
                for (uint32_t j = S->tao0[t] / TILE; j <= (S->tao0[t] + S->nc[t] - 1) / TILE; ++j) {
                    on[(size_t)MD_MIN(i, j) * S->nt + MD_MAX(i, j)] = 1;
                }
            }
        }
    }
    for (size_t i = 0; i < (size_t)S->nt * S->nt; ++i) c.tiles += on[i];
    c.ref_held = ref_pairs ? (double)ref_in / (double)ref_pairs : 1.0;
    free(on);
    return c;
}

// ---------------------------------------------------------------------------------------------------- D-weighted
// |phi_a(x)| <= f r^l sum_p |c_p| exp(-alpha_p r^2), r = |x - A|, f <= 1 (md_gto.h: f(i,j,k) = 1 / sqrt((2i-1)!!..)),
// |x^i y^j z^k| <= r^l. env[k] bounds it for every r >= k ENV_H: each term at max(r, its peak sqrt(l / (2 alpha))),
// non-increasing in k. For centres R apart, a point at distances a and b has a + b >= R; with env non-increasing,
// sup |phi_a phi_b| <= max over a in [0, R] of env_a(a) env_b(R - a), bounded cell by cell of the grid.

typedef struct dw_ctx_t {
    const shells_t* S;
    const double*   D;
    const double*   env;    // [ns][ENV_K]
    const uint64_t* row;    // [ns]: the first pair (s, s + 1) of row s in wmax and wsum
    float*          wmax;   // per pair s < t: max |D| over the pair's AOs times the bound on sup |phi_a phi_b|
    double*         wsum;   // per pair s < t: what the pair adds to rho anywhere, at most (both orders)
} dw_ctx_t;

static void dw_envelopes(double* env, const shells_t* S, const md_gto_basis_t* b) {
    for (uint32_t s = 0; s < S->ns; ++s) {
        const md_gto_shell_t* sh = &b->shells[s];
        const int l = (int)sh->l;
        for (int k = 0; k < ENV_K; ++k) {
            const double r = k * ENV_H;
            double v = 0.0;
            for (uint32_t p = 0; p < sh->num_primitives; ++p) {
                const double a = b->alpha[sh->primitive_offset + p];
                const double c = fabs((double)b->coeff[sh->primitive_offset + p]);
                const double rp = l > 0 ? sqrt((double)l / (2.0 * a)) : 0.0;
                const double x = r > rp ? r : rp;
                v += c * pow(x, l) * exp(-a * x * x);
            }
            env[(size_t)s * ENV_K + k] = v;
        }
    }
}

static inline double env_at(const double* e, double x) {
    if (x <= 0.0) return e[0];
    const double k = floor(x / ENV_H);
    return e[k >= ENV_K - 1 ? ENV_K - 1 : (int)k];
}

static double dw_product_bound(const double* ea, const double* eb, double R) {
    // cells [k h, (k + 1) h] of a: env_a(a) <= ea[k], env_b(R - a) <= env_b at R - (k + 1) h
    const int cells = (int)ceil(R / ENV_H);
    double m = 0.0;
    for (int k = 0; k < MD_MAX(cells, 1); ++k) {
        const double va = ea[k >= ENV_K - 1 ? ENV_K - 1 : k];
        const double vb = env_at(eb, R - (k + 1) * ENV_H);
        const double v = va * vb;
        if (v > m) m = v;
    }
    return m;
}

static void dw_rows(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    dw_ctx_t* c = (dw_ctx_t*)arg;
    const shells_t* S = c->S;
    const size_t N = S->nao;
    for (uint32_t s = beg; s < end; ++s) {
        uint64_t e = c->row[s];
        for (uint32_t t = s + 1; t < S->ns; ++t, ++e) {
            double dmax = 0.0, dsum = 0.0;
            for (uint32_t i = S->ao[s]; i < S->ao[s] + S->nc[s]; ++i) {
                for (uint32_t j = S->ao[t]; j < S->ao[t] + S->nc[t]; ++j) {
                    const double v = fabs(0.5 * (c->D[i * N + j] + c->D[j * N + i]));
                    dsum += v;
                    if (v > dmax) dmax = v;
                }
            }
            const double B = dw_product_bound(c->env + (size_t)s * ENV_K, c->env + (size_t)t * ENV_K, sqrt(shell_dist2(S, s, t)));
            c->wmax[e] = (float)(dmax * B);
            c->wsum[e] = 2.0 * dsum * B;
        }
    }
}

// ---------------------------------------------------------------------------------------------------- ESP

typedef struct esp_eval_t {
    const md_gto_int_charges_t* q;
    const float* pts;
    double* out;
} esp_eval_t;

static void esp_points(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    esp_eval_t* e = (esp_eval_t*)arg;
    md_gto_int_potential_xyz(e->out + beg, NULL, e->pts + 3 * (size_t)beg, end - beg, 0, e->q);
}

static double esp_eval(double* out, const md_gto_int_charges_t* q, const float* pts, uint32_t n, const md_executor_t* ex) {
    esp_eval_t e = { q, pts, out };
    const md_tick_t t0 = md_tick_now();
    md_parallel_for(ex, n, 1, esp_points, &e);
    return ms_since(t0);
}

// D (whole, or cut to a pattern: 0 outside it) as tiles that keep every pair, so that what is left out is D's zeros
static bool esp_build(md_gto_int_charges_t* out, const src_t* src, const double* D, double theta, double* ms) {
    md_gto_dmat_t dm;
    if (!qm_test_dmat_from_matrix(&dm, &src->basis, src->xyz, 0, D, md_get_heap_allocator())) return false;
    md_gto_int_charges_desc_t cd = {0};
    cd.basis = &src->basis;
    cd.atom_xyz = src->xyz;
    cd.dmat = &dm;
    cd.density_scale = -1.0;
    cd.point_xyz = src->xyz;
    cd.point_charge = src->Z;
    cd.num_points = src->num_atoms;
    cd.threshold = 1.0e-6;              // viamd's ESP_SCREENING_THRESHOLD
    cd.far_field_theta = theta;
    const md_tick_t t0 = md_tick_now();
    const bool ok = md_gto_int_charges_init(out, &cd, md_get_heap_allocator());
    *ms = ms_since(t0);
    md_gto_dmat_free(&dm);
    return ok;
}

// Points probe bohr from the nearest atom: around atoms taken in a stride, on a Fibonacci sphere each
static uint32_t esp_make_points(float* pts, uint32_t want, const src_t* src, double probe) {
    const uint32_t na = (uint32_t)src->num_atoms;
    const uint32_t dirs = 32;
    const double golden = 2.399963229728653;
    uint32_t n = 0;
    const uint32_t stride = 7919;   // a prime: the atoms visited are spread along the chain
    for (uint32_t k = 0; k < na * dirs && n < want; ++k) {
        const uint32_t a = (uint32_t)(((uint64_t)k * stride) % na);
        const uint32_t m = (k / na + k * 13) % dirs;
        const double z = 1.0 - (2.0 * m + 1.0) / dirs, r = sqrt(1.0 - z * z), ph = golden * m;
        double p[3] = { src->xyz[3 * a] + probe * r * cos(ph), src->xyz[3 * a + 1] + probe * r * sin(ph), src->xyz[3 * a + 2] + probe * z };
        double dmin2 = DBL_MAX;
        for (uint32_t b = 0; b < na; ++b) {
            const double dx = p[0] - src->xyz[3 * b], dy = p[1] - src->xyz[3 * b + 1], dz = p[2] - src->xyz[3 * b + 2];
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < dmin2) dmin2 = d2;
        }
        if (dmin2 < probe * probe * (1.0 - 1e-9)) continue;
        for (int c = 0; c < 3; ++c) pts[3 * n + c] = (float)p[c];
        n++;
    }
    return n;
}

static void mask_density(double* Dm, const double* D, const shells_t* S, const pattern_t* p) {
    const size_t N = S->nao;
    memcpy(Dm, D, sizeof(double) * N * N);
    for (uint32_t s = 0; s < S->ns; ++s) {
        for (uint32_t t = s + 1; t < S->ns; ++t) {
            if (pattern_has(p, s, t)) continue;
            for (uint32_t i = S->ao[s]; i < S->ao[s] + S->nc[s]; ++i) {
                for (uint32_t j = S->ao[t]; j < S->ao[t] + S->nc[t]; ++j) { Dm[i * N + j] = 0.0; Dm[j * N + i] = 0.0; }
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------------- main

typedef struct opts_t {
    double   deltas[MAX_DELTAS];
    int      num_deltas;
    uint32_t points, exact, threads;
    double   probe;
    bool     esp, topo;
    const char* data;
} opts_t;

static void run(const char* path, const opts_t* o, const md_executor_t* ex) {
    src_t src;
    printf("%s\n", path);
    if (!src_load(&src, path)) { printf("  could not load\n\n"); src_free(&src); return; }
    shells_t S;
    shells_make(&S, &src, 1e-6);    // one layout for every pattern: viamd's density dmat's
    const size_t N = src.num_ao;
    const double MB = 1.0 / (1024.0 * 1024.0);
    printf("  %zu atoms, %u shells, %zu AOs, %zu occupied orbitals\n", src.num_atoms, S.ns, N, src.num_occ);
    printf("  read: open %.0f ms, the dense matrix %.0f ms (%.1f MB), the orbitals %.0f ms (%.1f MB)\n",
           src.ms_open, src.ms_read_d, 8.0 * N * N * MB, src.ms_read_orbitals, 8.0 * src.num_occ * N * MB);
    const double tri = (double)N * (N + 1) / 2;
    const double tiles_total = (double)S.nt * (S.nt + 1) / 2;
    const double map_mb = 4.0 * S.nt * S.nt * MB;
    printf("  whole: dense double %.1f MB, triangle double %.1f MB, triangle float %.1f MB; %u x %u tiles of 16 (%.0f in the triangle), map %.2f MB\n\n",
           8.0 * N * N * MB, 8.0 * tri * MB, 4.0 * tri * MB, S.nt, S.nt, tiles_total, map_mb);

    // --- patterns
    enum { P_D6M0, P_D6M1, P_E7M1, P_GEOM };
    pattern_t geo[P_GEOM];
    double* r6 = calloc(S.ns, sizeof(double));
    double* r7 = calloc(S.ns, sizeof(double));
    for (uint32_t s = 0; s < S.ns; ++s) { r6[s] = shell_radius(&src.basis, s, 1e-6); r7[s] = shell_radius(&src.basis, s, 1e-7); }
    pattern_init(&geo[P_D6M0], "density 1e-6, margin 0", S.ns);   pattern_geometric(&geo[P_D6M0], &S, r6, NULL, 0.0);
    pattern_init(&geo[P_D6M1], "density 1e-6, margin 1 (viamd)", S.ns); pattern_geometric(&geo[P_D6M1], &S, r6, NULL, 1.0);
    pattern_init(&geo[P_E7M1], "ELF/LOL 1e-7, margin 1 (viamd)", S.ns); pattern_geometric(&geo[P_E7M1], &S, r7, NULL, 1.0);
    const pattern_t* ref = &geo[P_D6M1];

    pattern_t topo[4];
    int num_topo = 0;
    if (o->topo) {
        md_topo_gto_desc_t td = { .basis = &src.basis, .atom_xyz = src.xyz, .density_matrix = src.D, .executor = ex };
        md_topo_gto_reach_t reach;
        if (md_topo_gto_reach(&reach, &td, md_get_heap_allocator()) && reach.num_shells == S.ns) {
            printf("  topology (md_topo_gto_reach, setup %.0f ms): tau %.2e, rung %u, region margin %.2f bohr (CPU cubes alone %.2f)\n",
                   reach.ms_setup, reach.tau, reach.screening_rung, reach.region_margin, reach.region_margin_cpu);
            double rl_max = 0, rf_max = 0;
            for (uint32_t s = 0; s < S.ns; ++s) { rl_max = fmax(rl_max, reach.radius_local[s]); rf_max = fmax(rf_max, reach.radius_far[s]); }
            printf("    largest radius: local %.2f bohr, far %.2f bohr (the density's at 1e-6: ", rl_max, rf_max);
            double r6_max = 0; for (uint32_t s = 0; s < S.ns; ++s) r6_max = fmax(r6_max, r6[s]);
            printf("%.2f)\n\n", r6_max);
            pattern_init(&topo[0], "topology lists (local + region)", S.ns); pattern_geometric(&topo[0], &S, reach.radius_local, NULL, reach.region_margin);
            pattern_init(&topo[1], "topology polish (far)", S.ns);           pattern_geometric(&topo[1], &S, reach.radius_far, NULL, 0.0);
            pattern_init(&topo[2], "topology reads (both)", S.ns);
            for (size_t w = 0; w < (size_t)S.ns * topo[2].words; ++w) topo[2].bits[w] = topo[0].bits[w] | topo[1].bits[w];
            pattern_init(&topo[3], "topology GPU tiles today (far + region)", S.ns); pattern_geometric(&topo[3], &S, reach.radius_far, NULL, reach.region_margin);
            num_topo = 4;
            md_topo_gto_reach_free(&reach, md_get_heap_allocator());
        } else {
            printf("  topology: md_topo_gto_reach failed\n\n");
        }
    }

    // D-weighted: the bounds of every pair, then a pattern per delta
    const uint64_t num_pairs = (uint64_t)S.ns * (S.ns - 1) / 2;
    uint64_t* row = calloc(S.ns, sizeof(uint64_t));
    for (uint32_t s = 1; s < S.ns; ++s) row[s] = row[s - 1] + (S.ns - s);
    float*  wmax = calloc(MD_MAX(1, num_pairs), sizeof(float));
    double* wsum = calloc(MD_MAX(1, num_pairs), sizeof(double));
    pattern_t dw[MAX_DELTAS];
    {
        double* env = calloc((size_t)S.ns * ENV_K, sizeof(double));
        dw_envelopes(env, &S, &src.basis);
        dw_ctx_t c = { &S, src.D, env, row, wmax, wsum };
        const md_tick_t t0 = md_tick_now();
        md_parallel_for(ex, S.ns, 8, dw_rows, &c);
        const double ms = ms_since(t0);
        free(env);
        for (int d = 0; d < o->num_deltas; ++d) {
            char name[96];
            snprintf(name, sizeof(name), "D-weighted, delta %.0e", o->deltas[d]);
            pattern_init(&dw[d], name, S.ns);
            for (uint32_t s = 0; s < S.ns; ++s) {
                pattern_set(&dw[d], s, s);
                for (uint32_t t = s + 1; t < S.ns; ++t) if (wmax[row[s] + (t - s - 1)] > o->deltas[d]) pattern_set(&dw[d], s, t);
            }
        }
        double E_all = 0.0;
        for (uint64_t e = 0; e < num_pairs; ++e) E_all += wsum[e];
        printf("  D-weighted bounds: %.0f ms (envelopes on %.2f bohr to %.0f bohr; rho only); every pair between different shells together: %.2e\n\n",
               ms, ENV_H, ENV_H * (ENV_K - 1), E_all);
    }

    // --- the census
    printf("  %-42s %10s %12s %7s %9s %9s %9s %10s %11s %11s\n", "pattern", "shell prs", "AO pairs", "% tri", "tiles", "% tiles", "MB float", "MB double", "viamd's held", "E left out");
    const pattern_t* all[P_GEOM + 4 + MAX_DELTAS];
    int na = 0;
    for (int i = 0; i < P_GEOM; ++i) all[na++] = &geo[i];
    for (int i = 0; i < num_topo; ++i) all[na++] = &topo[i];
    for (int i = 0; i < o->num_deltas; ++i) all[na++] = &dw[i];
    for (int i = 0; i < na; ++i) {
        const census_t c = pattern_census(all[i], &S, ref, row, wsum);
        printf("  %-42s %10llu %12llu %6.1f%% %9llu %8.1f%% %9.1f %10.1f %10.2f%%",
               all[i]->name, (unsigned long long)c.shell_pairs, (unsigned long long)c.ao_pairs, 100.0 * c.ao_pairs / tri,
               (unsigned long long)c.tiles, 100.0 * c.tiles / tiles_total, 1024.0 * c.tiles * MB + map_mb, 2048.0 * c.tiles * MB + map_mb, 100.0 * c.ref_held);
        printf("  %9.1e\n", c.E);
        fflush(stdout);
    }
    printf("  (MB: the tiles and the %.2f MB map. E: a bound on |rho - rho'| anywhere from the shell pairs left out, e/bohr^3,\n"
           "   their bounds summed: pessimistic for a large system, where most of them are nowhere near one another)\n\n", map_mb);

    // --- ESP
    if (o->esp) {
        float* pts = calloc(3 * (size_t)MD_MAX(1, o->points), sizeof(float));
        const uint32_t np = esp_make_points(pts, o->points, &src, o->probe);
        md_gto_int_charges_t full;
        double ms_full = 0.0;
        if (np == 0 || !esp_build(&full, &src, src.D, MD_GTO_INT_FAR_FIELD_THETA, &ms_full)) {
            printf("  ESP: %s\n\n", np == 0 ? "no points" : "the charge distribution failed");
        } else {
            double* V0 = calloc(np, sizeof(double));
            double* V = calloc(np, sizeof(double));
            const double ms_eval = esp_eval(V0, &full, pts, np, ex);
            const md_gto_int_moments_t m0 = md_gto_int_charges_moments(&full, NULL);
            double rms0 = 0.0, max0 = 0.0;
            for (uint32_t i = 0; i < np; ++i) { rms0 += V0[i] * V0[i]; max0 = fmax(max0, fabs(V0[i])); }
            rms0 = sqrt(rms0 / np);
            printf("  ESP (threshold 1e-6, far field %.2f), %u points %.1f bohr from the nearest atom: the whole D in %.0f ms, %u gaussians, net charge %+.4f e;\n"
                   "  |V| at the points: RMS %.3e, max %.3e au (%.0f ms for the points)\n",
                   MD_GTO_INT_FAR_FIELD_THETA, np, o->probe, ms_full, full.num_gaussians, m0.charge, rms0, max0, ms_eval);
            if (o->exact > 0) {
                md_gto_int_charges_t ex_q;
                double ms_x = 0.0;
                const uint32_t nx = MD_MIN(o->exact, np);
                if (esp_build(&ex_q, &src, src.D, 0.0, &ms_x)) {
                    const double ms_e = esp_eval(V, &ex_q, pts, nx, ex);
                    double mx = 0.0;
                    for (uint32_t i = 0; i < nx; ++i) mx = fmax(mx, fabs(V[i] - V0[i]));
                    printf("  the far field's own error at %u points: max %.2e au (without it %.0f ms for them)\n", nx, mx, ms_e);
                    md_gto_int_charges_free(&ex_q, md_get_heap_allocator());
                }
            }
            printf("\n  %-42s %10s %9s %11s %11s %11s\n", "D cut to", "gaussians", "build ms", "dQ (e)", "max |dV|", "RMS dV");
            double* Dm = md_alloc(md_get_heap_allocator(), sizeof(double) * N * N);
            const pattern_t* esp_pats[2 + MAX_DELTAS];
            int ne = 0;
            esp_pats[ne++] = &geo[P_D6M1];
            esp_pats[ne++] = &geo[P_E7M1];
            for (int i = 0; i < o->num_deltas; ++i) esp_pats[ne++] = &dw[i];
            for (int i = 0; i < ne; ++i) {
                mask_density(Dm, src.D, &S, esp_pats[i]);
                md_gto_int_charges_t q;
                double ms = 0.0;
                if (!esp_build(&q, &src, Dm, MD_GTO_INT_FAR_FIELD_THETA, &ms)) { printf("  %-42s failed\n", esp_pats[i]->name); continue; }
                esp_eval(V, &q, pts, np, ex);
                const md_gto_int_moments_t m = md_gto_int_charges_moments(&q, NULL);
                double mx = 0.0, rms = 0.0;
                for (uint32_t k = 0; k < np; ++k) { const double dv = V[k] - V0[k]; mx = fmax(mx, fabs(dv)); rms += dv * dv; }
                printf("  %-42s %10u %9.0f %+11.2e %11.2e %11.2e\n", esp_pats[i]->name, q.num_gaussians, ms, m.charge - m0.charge, mx, sqrt(rms / np));
                fflush(stdout);
                md_gto_int_charges_free(&q, md_get_heap_allocator());
            }
            md_free(md_get_heap_allocator(), Dm, sizeof(double) * N * N);
            free(V0); free(V);
            md_gto_int_charges_free(&full, md_get_heap_allocator());
            printf("\n");
        }
        free(pts);
    }

    for (int i = 0; i < P_GEOM; ++i) pattern_free(&geo[i]);
    for (int i = 0; i < num_topo; ++i) pattern_free(&topo[i]);
    for (int i = 0; i < o->num_deltas; ++i) pattern_free(&dw[i]);
    free(r6); free(r7);
    free(row); free(wmax); free(wsum);
    shells_free(&S);
    src_free(&src);
}

static void quiet_log(struct md_logger_o* inst, md_log_type_t type, const char* msg) {
    (void)inst;
    if (type == MD_LOG_TYPE_ERROR) fprintf(stderr, "%s\n", msg);
}
static md_logger_i quiet_logger = { NULL, quiet_log };

static bool file_exists(const char* p) {
    FILE* f = fopen(p, "rb");
    if (f) fclose(f);
    return f != NULL;
}

int main(int argc, char** argv) {
    opts_t o = { .points = 256, .exact = 16, .probe = 4.0, .esp = true, .topo = true, .data = MD_BENCHMARK_DATA_DIR };
    const double def_deltas[] = { 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16 };
    o.num_deltas = 6;
    memcpy(o.deltas, def_deltas, sizeof(def_deltas));
    const char* files[32];
    int num_files = 0;
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        const bool more = i + 1 < argc;
        if (!strcmp(a, "--deltas") && more) {
            o.num_deltas = 0;
            char buf[512];
            snprintf(buf, sizeof(buf), "%s", argv[++i]);
            for (char* t = strtok(buf, ","); t && o.num_deltas < MAX_DELTAS; t = strtok(NULL, ",")) o.deltas[o.num_deltas++] = atof(t);
        }
        else if (!strcmp(a, "--points") && more)  o.points = (uint32_t)atoi(argv[++i]);
        else if (!strcmp(a, "--exact") && more)   o.exact = (uint32_t)atoi(argv[++i]);
        else if (!strcmp(a, "--probe") && more)   o.probe = atof(argv[++i]);
        else if (!strcmp(a, "--threads") && more) o.threads = (uint32_t)atoi(argv[++i]);
        else if (!strcmp(a, "--data") && more)    o.data = argv[++i];
        else if (!strcmp(a, "--no-esp"))         o.esp = false;
        else if (!strcmp(a, "--no-topo"))        o.topo = false;
        else if (a[0] == '-') {
            printf("usage: %s [--deltas a,b,..] [--points n] [--exact n] [--probe bohr] [--no-esp] [--no-topo] [--threads n] [--data dir] [file ...]\n", argv[0]);
            return 1;
        }
        else if (num_files < 32) files[num_files++] = a;
    }
    if (default_logger) {
        md_log_unregister(default_logger);
        md_log_register(&quiet_logger);
    }

    md_os_sys_info_t si = {0};
    md_os_sys_info_query(&si);
    test_executor_t exec = test_executor_threads(o.threads > 0 ? o.threads : (si.num_virtual_cores > 0 ? si.num_virtual_cores : 1));
    const md_executor_t* ex = test_executor(&exec);

    printf("%s | source %s | %u threads\n", MD_DMAT_CENSUS_BENCH_REV, MD_BENCH_SOURCE_REV, ex->num_threads);
    printf("deltas:");
    for (int d = 0; d < o.num_deltas; ++d) printf(" %.0e", o.deltas[d]);
    printf("\n\n");

    char path[1024];
    if (num_files == 0) {
        snprintf(path, sizeof(path), "%s/vlx/protein.h5", o.data);
        if (!file_exists(path)) snprintf(path, sizeof(path), "%s/vlx/mol.h5", o.data);
        run(path, &o, ex);
    }
    for (int i = 0; i < num_files; ++i) {
        if (file_exists(files[i])) snprintf(path, sizeof(path), "%s", files[i]);
        else snprintf(path, sizeof(path), "%s/vlx/%s", o.data, files[i]);
        run(path, &o, ex);
    }
    return 0;
}

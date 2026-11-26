// md_basis_census_bench
//
// What the density's GPU forms would cost for a molecule in other basis sets, counted on the CPU (no GPU). The tile
// pattern of md_gto_gpu_dmat depends only on the geometry, the basis and the cutoff, never on the coefficients, and
// so do the AOs each block of the GEMM path screens in: a file needs only a geometry and a basis set name, no SCF.
// (benchmark/basis_census_inputs.py writes such files, with the basis sets, from the geometry of an existing one.)
//
//   AOs       Cartesian (what mdlib evaluates) and spherical (what VeloxChem counts), and the occupied orbitals
//             C_occ in float in each, n_occ from the file's alpha electrons (a closed shell)
//   blocks    the AOs that reach each block of 8 x 8 x 8 voxels (the GEMM path's: AABB within the shell's radius),
//             over the GEMM benchmark's grid (the atoms' box +- --margin, --h spacing), averaged over the non-empty
//             ones, and the GEMM's multiply-adds per voxel in the D form (n (n + 1) / 2) and in the orbital form
//             (n n_occ), in both bases
//   tiles     the 16 x 16 tile pairs (I <= J) a pattern stores, the fraction of all, and their bytes in float, for
//             the density's pattern (cutoff 1e-6, margin 1 bohr: viamd's) and ELF/LOL's (1e-7, margin 1), in six
//             AO orders, Cartesian and spherical. The tiles beat C_occ in memory below the fraction 2 n_occ / N.
//               file      VeloxChem's matrix rows: angular momentum, component, atom, function
//               shells    mdlib's shell order: angular momentum, atom, function, component
//               tiled     md_gto_gpu_dmat's today: radius class (2 bohr), atom (the file's order), shell
//               atom      atom (the file's order), then its shells by radius, longest first
//               hilbert   atoms along a Hilbert curve, then each atom's shells by radius, longest first
//               cls+hilb  radius class (2 bohr), then atoms along the Hilbert curve, then shell
//
//   md_basis_census_bench [--h bohr] [--margin bohr] [--orders] file.h5 ...
//
//   --orders  every order (default: tiled only, the others for the density's pattern alone)

#include <md_vlx.h>
#include <md_system.h>
#include <md_gto.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_unit.h>

#include "bench_rev.h"      // MD_BENCH_SOURCE_REV, generated at build time by bench_version.cmake

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Bumped with every change to this program's measurement or output
#define MD_BASIS_CENSUS_BENCH_REV "basis-census r1"

#define ANGSTROM_TO_BOHR 1.8897261246257702
#define TILE 16u
#define BLOCK 8u

static double ms_since(md_tick_t t0) { return md_tick_to_milliseconds(md_tick_now() - t0); }

// ---------------------------------------------------------------------------------------------------- the input

typedef struct src_t {
    md_allocator_i* arena;
    md_gto_basis_t  basis;
    char            basis_name[64];
    float*          xyz;            // bohr
    uint32_t        num_atoms, ns;
    uint32_t        n_occ;
    uint32_t*       atom_first;     // the shells of atom a: atom_shell[atom_first[a] .. atom_first[a + 1])
    uint32_t*       atom_shell;
} src_t;

static bool src_load(src_t* s, const char* path) {
    memset(s, 0, sizeof(*s));
    s->arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(64));
    md_system_t sys = { .alloc = s->arena };
    md_system_state_t state = { .alloc = s->arena };
    if (!md_vlx_system_init_from_file(&sys, &state, md_str_from_cstr(path), NULL)) return false;
    const md_attribute_t* c = md_attributes_find(&sys.attributes, MD_STR_LIT("qm/atom/coordinate"));
    const md_attribute_t* e = md_attributes_find(&sys.attributes, MD_STR_LIT("vlx/electron_count/alpha"));
    const md_attribute_t* b = md_attributes_find(&sys.attributes, MD_STR_LIT("vlx/basis_set"));
    if (!c || !e || !md_gto_basis_extract_attributes(&s->basis, &sys.attributes, s->arena)) {
        fprintf(stderr, "%s: no basis, geometry or electron count\n", path);
        return false;
    }
    if (b) {
        const md_str_t name = md_attribute_str(&sys.attributes, b, 0);
        snprintf(s->basis_name, sizeof(s->basis_name), MD_STR_FMT, MD_STR_ARG(name));
    }
    s->num_atoms = (uint32_t)md_attribute_value_count(&c->format);
    double* xyz = md_arena_allocator_push(s->arena, sizeof(double) * 3 * s->num_atoms);
    md_attribute_extract_f64(xyz, 3 * s->num_atoms, &sys.attributes, c, md_attribute_slice_all(), md_unit_none());
    s->xyz = md_arena_allocator_push(s->arena, sizeof(float) * 3 * s->num_atoms);
    for (size_t i = 0; i < 3 * (size_t)s->num_atoms; ++i) s->xyz[i] = (float)(xyz[i] * ANGSTROM_TO_BOHR);
    double ne = 0.0;
    md_attribute_extract_f64(&ne, 1, &sys.attributes, e, md_attribute_slice_all(), md_unit_none());
    s->n_occ = (uint32_t)ne;
    s->ns = s->basis.num_shells;
    s->atom_first = md_arena_allocator_push(s->arena, sizeof(uint32_t) * (s->num_atoms + 1));
    s->atom_shell = md_arena_allocator_push(s->arena, sizeof(uint32_t) * MD_MAX(1, s->ns));
    memset(s->atom_first, 0, sizeof(uint32_t) * (s->num_atoms + 1));
    for (uint32_t i = 0; i < s->ns; ++i) s->atom_first[s->basis.shells[i].atom_idx + 1]++;
    for (uint32_t a = 0; a < s->num_atoms; ++a) s->atom_first[a + 1] += s->atom_first[a];
    uint32_t* fill = calloc(s->num_atoms, sizeof(uint32_t));
    for (uint32_t i = 0; i < s->ns; ++i) {
        const uint32_t a = s->basis.shells[i].atom_idx;
        s->atom_shell[s->atom_first[a] + fill[a]++] = i;
    }
    free(fill);
    return true;
}

static void src_free(src_t* s) {
    if (s->arena) md_arena_allocator_destroy(s->arena);
    memset(s, 0, sizeof(*s));
}

static inline uint32_t num_comp(uint32_t l, bool sph) { return sph ? md_gto_num_sph_ao(l) : md_gto_num_cart_ao(l); }

// A shell's screening radius at a cutoff, as md_gto_gpu_basis and md_gto_gpu_dmat make it (gto_expand_basis_gpu_shells:
// the largest over components and primitives, the coefficient rounded to float); as bench_dmat_census.c
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

// ---------------------------------------------------------------------------------------------------- AO orders

enum { ORDER_FILE, ORDER_SHELLS, ORDER_TILED, ORDER_ATOM, ORDER_HILBERT, ORDER_CLS_HILB, ORDER_COUNT };
static const char* order_name[ORDER_COUNT] = { "file", "shells", "tiled", "atom", "hilbert", "cls+hilb" };

// The index along a 3D Hilbert curve of a point on a 2^bits grid (Skilling, "Programming the Hilbert curve", 2004)
static uint64_t hilbert3(uint32_t x, uint32_t y, uint32_t z, int bits) {
    uint32_t X[3] = { x, y, z };
    const uint32_t M = 1u << (bits - 1);
    for (uint32_t Q = M; Q > 1; Q >>= 1) {
        const uint32_t P = Q - 1;
        for (int i = 0; i < 3; ++i) {
            if (X[i] & Q) X[0] ^= P;
            else { const uint32_t t = (X[0] ^ X[i]) & P; X[0] ^= t; X[i] ^= t; }
        }
    }
    for (int i = 1; i < 3; ++i) X[i] ^= X[i - 1];
    uint32_t t = 0;
    for (uint32_t Q = M; Q > 1; Q >>= 1) if (X[2] & Q) t ^= Q - 1;
    for (int i = 0; i < 3; ++i) X[i] ^= t;
    uint64_t h = 0;
    for (int b = bits - 1; b >= 0; --b) for (int i = 0; i < 3; ++i) h = (h << 1) | ((X[i] >> b) & 1u);
    return h;
}

static int cmp_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return x < y ? -1 : x > y;
}

// Each atom's rank along the Hilbert curve through the atoms' box
static void hilbert_ranks(uint32_t* rank, const src_t* s) {
    float lo[3] = { 1e30f, 1e30f, 1e30f }, hi[3] = { -1e30f, -1e30f, -1e30f };
    for (uint32_t a = 0; a < s->num_atoms; ++a) for (int k = 0; k < 3; ++k) {
        lo[k] = MD_MIN(lo[k], s->xyz[3 * a + k]);
        hi[k] = MD_MAX(hi[k], s->xyz[3 * a + k]);
    }
    const float ext = MD_MAX(hi[0] - lo[0], MD_MAX(hi[1] - lo[1], hi[2] - lo[2])) + 1e-3f;
    uint64_t* key = malloc(sizeof(uint64_t) * MD_MAX(1, s->num_atoms));
    for (uint32_t a = 0; a < s->num_atoms; ++a) {
        uint32_t q[3];
        for (int k = 0; k < 3; ++k) q[k] = (uint32_t)((s->xyz[3 * a + k] - lo[k]) / ext * 1023.0f);
        key[a] = (hilbert3(q[0], q[1], q[2], 10) << 20) | a;
    }
    qsort(key, s->num_atoms, sizeof(uint64_t), cmp_u64);
    for (uint32_t i = 0; i < s->num_atoms; ++i) rank[key[i] & 0xFFFFFu] = i;
    free(key);
}

typedef struct kv_t { uint64_t k; uint32_t i; } kv_t;

static int cmp_kv(const void* x, const void* y) {
    const uint64_t a = ((const kv_t*)x)->k, b = ((const kv_t*)y)->k;
    return a < b ? -1 : a > b;
}

// Per shell, the distinct tiles its AOs fall in under an order: tiles[first[s] .. first[s + 1])
typedef struct shell_tiles_t {
    uint32_t  nao, nt;
    uint32_t* first;
    uint32_t* tiles;
} shell_tiles_t;

static void shell_tiles_free(shell_tiles_t* t) { free(t->first); free(t->tiles); memset(t, 0, sizeof(*t)); }

// r: the radii the order's classes are taken from (the pattern's own, as md_gto_gpu_dmat_create does)
static void order_make(shell_tiles_t* out, const src_t* s, int order, bool sph, const float* r, const uint32_t* hrank) {
    const md_gto_basis_t* b = &s->basis;
    const uint32_t ns = s->ns;
    uint32_t* ao0 = malloc(sizeof(uint32_t) * (ns + 1));
    uint32_t* func = malloc(sizeof(uint32_t) * MD_MAX(1, ns));     // the shell's rank among its atom's of its l
    ao0[0] = 0;
    for (uint32_t i = 0; i < ns; ++i) ao0[i + 1] = ao0[i] + num_comp(b->shells[i].l, sph);
    const uint32_t nao = ao0[ns];
    {
        // shells come angular momentum first, then atom, then function (md_gto.h): the rank within (l, atom) is a count
        uint32_t prev_l = UINT32_MAX, prev_a = UINT32_MAX, k = 0;
        for (uint32_t i = 0; i < ns; ++i) {
            if (b->shells[i].l != prev_l || b->shells[i].atom_idx != prev_a) k = 0;
            prev_l = b->shells[i].l; prev_a = b->shells[i].atom_idx;
            func[i] = k++;
        }
    }
    uint64_t* key = malloc(sizeof(uint64_t) * MD_MAX(1, nao));
    for (uint32_t i = 0; i < ns; ++i) {
        const uint64_t l = b->shells[i].l, a = b->shells[i].atom_idx;
        const uint64_t cls = (uint64_t)(r[i] * 0.5f);
        const uint64_t inv = 0xFFFFu - MD_MIN(0xFFFFu, (uint32_t)(r[i] * 256.0f));
        for (uint32_t c = 0; c < ao0[i + 1] - ao0[i]; ++c) {
            uint64_t k = 0;
            switch (order) {
            case ORDER_FILE:     k = (l << 40) | ((uint64_t)c << 32) | (a << 8) | func[i]; break;
            case ORDER_SHELLS:   k = ((uint64_t)i << 4) | c; break;
            case ORDER_TILED:    k = (cls << 44) | (a << 24) | ((uint64_t)i << 4) | c; break;
            case ORDER_ATOM:     k = (a << 40) | (inv << 24) | ((uint64_t)i << 4) | c; break;
            case ORDER_HILBERT:  k = ((uint64_t)hrank[a] << 40) | (inv << 24) | ((uint64_t)i << 4) | c; break;
            case ORDER_CLS_HILB: k = (cls << 44) | ((uint64_t)hrank[a] << 24) | ((uint64_t)i << 4) | c; break;
            }
            key[ao0[i] + c] = k;
        }
    }
    // the AOs sorted by key (the keys are unique: each holds its shell or its atom, l, component and function)
    kv_t* kv = malloc(sizeof(kv_t) * MD_MAX(1, nao));
    for (uint32_t i = 0; i < nao; ++i) kv[i] = (kv_t){ key[i], i };
    qsort(kv, nao, sizeof(kv_t), cmp_kv);
    uint32_t* pos = malloc(sizeof(uint32_t) * MD_MAX(1, nao));
    for (uint32_t p = 0; p < nao; ++p) pos[kv[p].i] = p;
    out->nao = nao;
    out->nt = (nao + TILE - 1) / TILE;
    out->first = malloc(sizeof(uint32_t) * (ns + 1));
    out->tiles = malloc(sizeof(uint32_t) * MD_MAX(1, nao));
    uint32_t n = 0;
    for (uint32_t i = 0; i < ns; ++i) {
        out->first[i] = n;
        for (uint32_t ao = ao0[i]; ao < ao0[i + 1]; ++ao) {
            const uint32_t t = pos[ao] / TILE;
            bool seen = false;
            for (uint32_t j = out->first[i]; j < n; ++j) seen |= out->tiles[j] == t;
            if (!seen) out->tiles[n++] = t;
        }
    }
    out->first[ns] = n;
    free(ao0); free(func); free(key); free(kv); free(pos);
}

// ---------------------------------------------------------------------------------------------------- patterns

// The tile pairs I <= J of the shell pairs whose centres are no further apart than their radii plus margin (and both
// radii > 0), as md_gto_gpu_dmat_create; visited atom pair by atom pair
static uint64_t pattern_tiles(const src_t* s, const shell_tiles_t* st, const float* r, double margin, uint8_t* on) {
    const uint32_t nt = st->nt;
    memset(on, 0, (size_t)nt * nt);
    float* ra = calloc(s->num_atoms, sizeof(float));
    for (uint32_t i = 0; i < s->ns; ++i) {
        const uint32_t a = s->basis.shells[i].atom_idx;
        ra[a] = MD_MAX(ra[a], r[i]);
    }
    for (uint32_t a = 0; a < s->num_atoms; ++a) {
        const float* p = s->xyz + 3 * a;
        for (uint32_t c = a; c < s->num_atoms; ++c) {
            const float* q = s->xyz + 3 * c;
            const double dx = p[0] - q[0], dy = p[1] - q[1], dz = p[2] - q[2];
            const double d2 = dx * dx + dy * dy + dz * dz;
            const double am = (double)ra[a] + (double)ra[c] + margin;
            if (d2 > am * am) continue;
            for (uint32_t x = s->atom_first[a]; x < s->atom_first[a + 1]; ++x) {
                const uint32_t si = s->atom_shell[x];
                if (!(r[si] > 0.0f)) continue;
                for (uint32_t y = (a == c ? x : s->atom_first[c]); y < s->atom_first[c + 1]; ++y) {
                    const uint32_t ui = s->atom_shell[y];
                    if (!(r[ui] > 0.0f)) continue;
                    const double lim = (double)r[si] + (double)r[ui] + margin;
                    if (d2 > lim * lim) continue;
                    for (uint32_t i = st->first[si]; i < st->first[si + 1]; ++i) {
                        for (uint32_t j = st->first[ui]; j < st->first[ui + 1]; ++j) {
                            const uint32_t I = st->tiles[i], J = st->tiles[j];
                            on[(size_t)MD_MIN(I, J) * nt + MD_MAX(I, J)] = 1;
                        }
                    }
                }
            }
        }
    }
    free(ra);
    uint64_t n = 0;
    for (size_t i = 0; i < (size_t)nt * nt; ++i) n += on[i];
    return n;
}

// ---------------------------------------------------------------------------------------------------- blocks

typedef struct block_stats_t {
    uint64_t blocks, nonempty;
    double   sum_n[2], sum_d[2], sum_c[2], max_n[2];   // [Cartesian, spherical]
} block_stats_t;

static inline double aabb_dist2(const float* p, const float* lo, const float* hi) {
    double d2 = 0.0;
    for (int k = 0; k < 3; ++k) {
        const double v = p[k] < lo[k] ? lo[k] - p[k] : (p[k] > hi[k] ? p[k] - hi[k] : 0.0);
        d2 += v * v;
    }
    return d2;
}

// The AOs that reach each block of BLOCK^3 voxels of the GEMM benchmark's grid (eval_gto_density_gemm.slang's
// shell_passes: the block's AABB within the shell's radius)
static block_stats_t block_census(const src_t* s, const float* r, double h, double margin) {
    block_stats_t st = {0};
    float lo[3] = { 1e30f, 1e30f, 1e30f }, hi[3] = { -1e30f, -1e30f, -1e30f };
    for (uint32_t a = 0; a < s->num_atoms; ++a) for (int k = 0; k < 3; ++k) {
        lo[k] = MD_MIN(lo[k], s->xyz[3 * a + k]);
        hi[k] = MD_MAX(hi[k], s->xyz[3 * a + k]);
    }
    uint32_t nb[3];
    for (int k = 0; k < 3; ++k) {
        lo[k] -= (float)margin; hi[k] += (float)margin;
        const uint32_t dim = (uint32_t)ceil((hi[k] - lo[k]) / h);
        nb[k] = (dim + BLOCK - 1) / BLOCK;
    }
    float* ra = calloc(s->num_atoms, sizeof(float));
    for (uint32_t i = 0; i < s->ns; ++i) {
        const uint32_t a = s->basis.shells[i].atom_idx;
        ra[a] = MD_MAX(ra[a], r[i]);
    }
    const float edge = (float)(h * BLOCK);
    for (uint32_t z = 0; z < nb[2]; ++z) for (uint32_t y = 0; y < nb[1]; ++y) for (uint32_t x = 0; x < nb[0]; ++x) {
        const float bl[3] = { lo[0] + x * edge, lo[1] + y * edge, lo[2] + z * edge };
        const float bh[3] = { bl[0] + edge, bl[1] + edge, bl[2] + edge };
        uint64_t n[2] = { 0, 0 };
        for (uint32_t a = 0; a < s->num_atoms; ++a) {
            const double d2 = aabb_dist2(s->xyz + 3 * a, bl, bh);
            if (d2 >= (double)ra[a] * ra[a]) continue;
            for (uint32_t i = s->atom_first[a]; i < s->atom_first[a + 1]; ++i) {
                const uint32_t si = s->atom_shell[i];
                if (r[si] > 0.0f && d2 < (double)r[si] * r[si]) {
                    n[0] += num_comp(s->basis.shells[si].l, false);
                    n[1] += num_comp(s->basis.shells[si].l, true);
                }
            }
        }
        st.blocks++;
        if (!n[0]) continue;
        st.nonempty++;
        for (int k = 0; k < 2; ++k) {
            st.sum_n[k] += (double)n[k];
            st.sum_d[k] += (double)n[k] * (double)(n[k] + 1) * 0.5;
            st.sum_c[k] += (double)n[k] * (double)s->n_occ;
            st.max_n[k] = MD_MAX(st.max_n[k], (double)n[k]);
        }
    }
    free(ra);
    return st;
}

// ---------------------------------------------------------------------------------------------------- main

static double mb(double bytes) { return bytes / (1024.0 * 1024.0); }

int main(int argc, char** argv) {
    double h = 0.25, margin = 4.0;
    bool all_orders = false;
    const char* files[64];
    int num_files = 0;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--h") && i + 1 < argc) h = atof(argv[++i]);
        else if (!strcmp(argv[i], "--margin") && i + 1 < argc) margin = atof(argv[++i]);
        else if (!strcmp(argv[i], "--orders")) all_orders = true;
        else if (argv[i][0] == '-') {
            printf("usage: %s [--h bohr] [--margin bohr] [--orders] file.h5 ...\n", argv[0]);
            return 1;
        } else if (num_files < 64) files[num_files++] = argv[i];
    }
    printf("md_basis_census_bench %s, mdlib %s\n", MD_BASIS_CENSUS_BENCH_REV, MD_BENCH_SOURCE_REV);
    printf("blocks: %u^3 voxels of %.3g bohr over the atoms' box +- %.3g bohr, the density's cutoff 1e-6 | tiles: %ux%u, float\n\n",
           BLOCK, h, margin, TILE, TILE);

    static const struct { const char* name; double cutoff, margin; } pats[] = {
        { "density 1e-6 +1", 1.0e-6, 1.0 },
        { "ELF/LOL 1e-7 +1", 1.0e-7, 1.0 },
    };
    for (int fi = 0; fi < num_files; ++fi) {
        src_t s;
        if (!src_load(&s, files[fi])) { src_free(&s); continue; }
        const md_tick_t t0 = md_tick_now();
        size_t nao[2] = { md_gto_basis_num_ao(&s.basis), md_gto_basis_num_sph_ao(&s.basis) };
        uint32_t max_l = 0;
        for (uint32_t i = 0; i < s.ns; ++i) max_l = MD_MAX(max_l, s.basis.shells[i].l);
        printf("%s: %s, %u atoms, %u shells (l <= %u), n_occ %u\n", files[fi], s.basis_name, s.num_atoms, s.ns, max_l, s.n_occ);
        printf("  AOs: %zu Cartesian, %zu spherical (%.3fx) | AOs per atom %.1f / %.1f\n", nao[0], nao[1], (double)nao[0] / (double)nao[1],
               (double)nao[0] / s.num_atoms, (double)nao[1] / s.num_atoms);
        printf("  C_occ in float: %.1f MB Cartesian, %.1f MB spherical | the tiles beat it below the fraction %.4f / %.4f\n",
               mb(4.0 * s.n_occ * nao[0]), mb(4.0 * s.n_occ * nao[1]), 2.0 * s.n_occ / nao[0], 2.0 * s.n_occ / nao[1]);
        printf("  the whole upper triangle in float: %.1f MB Cartesian, %.1f MB spherical\n",
               mb(2.0 * nao[0] * (nao[0] + 1.0)), mb(2.0 * nao[1] * (nao[1] + 1.0)));

        float* r = malloc(sizeof(float) * MD_MAX(1, s.ns));
        uint32_t* hrank = malloc(sizeof(uint32_t) * MD_MAX(1, s.num_atoms));
        hilbert_ranks(hrank, &s);

        for (uint32_t i = 0; i < s.ns; ++i) r[i] = shell_radius(&s.basis, i, 1.0e-6);
        const block_stats_t bs = block_census(&s, r, h, margin);
        for (int k = 0; k < 2; ++k) {
            const double ne = (double)MD_MAX(1, bs.nonempty);
            printf("  blocks (%s): %llu of %llu non-empty, AOs per non-empty block %.0f (max %.0f) | multiply-adds per voxel: "
                   "D form %.0f, orbital form %.0f (D / orbital %.3f)\n",
                   k ? "spherical" : "Cartesian", (unsigned long long)bs.nonempty, (unsigned long long)bs.blocks, bs.sum_n[k] / ne,
                   bs.max_n[k], bs.sum_d[k] / ne, bs.sum_c[k] / ne, bs.sum_c[k] > 0 ? bs.sum_d[k] / bs.sum_c[k] : 0.0);
        }

        printf("  %-16s %-9s | %10s %8s %10s | %10s %8s %10s\n", "pattern", "order", "Cart tiles", "fraction", "MB", "sph tiles",
               "fraction", "MB");
        for (size_t p = 0; p < MD_COUNTOF(pats); ++p) {
            for (uint32_t i = 0; i < s.ns; ++i) r[i] = shell_radius(&s.basis, i, pats[p].cutoff);
            for (int o = 0; o < ORDER_COUNT; ++o) {
                if (!all_orders && p > 0 && o != ORDER_TILED) continue;
                uint64_t tiles[2];
                double frac[2];
                for (int k = 0; k < 2; ++k) {
                    shell_tiles_t st;
                    order_make(&st, &s, o, k == 1, r, hrank);
                    uint8_t* on = malloc((size_t)st.nt * st.nt);
                    tiles[k] = pattern_tiles(&s, &st, r, pats[p].margin, on);
                    frac[k] = (double)tiles[k] / (0.5 * st.nt * (st.nt + 1.0));
                    free(on);
                    shell_tiles_free(&st);
                }
                printf("  %-16s %-9s | %10llu %8.4f %10.1f | %10llu %8.4f %10.1f\n", o == 0 || (!all_orders && p > 0) ? pats[p].name : "",
                       order_name[o], (unsigned long long)tiles[0], frac[0], mb(1024.0 * tiles[0]), (unsigned long long)tiles[1], frac[1],
                       mb(1024.0 * tiles[1]));
            }
        }
        printf("  (%.1f s)\n\n", ms_since(t0) * 1e-3);
        free(r); free(hrank);
        src_free(&s);
    }
    return 0;
}

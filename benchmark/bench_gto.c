#include "ubench.h"

#include <md_vlx.h>
#include <md_system.h>
#include <md_gto.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#endif

#include <float.h>
#include <math.h>
#include <stdio.h>

// Molecular orbital (the HOMO) and total electron density of test_data/vlx/mol.h5 (26 atoms, def2-SVP,
// B3LYP ground state), on a 0.4 bohr grid reaching 6 bohr past the atoms - the same dataset and grid as
// bench_gto_int.
//
// GPU benchmarks evaluate the whole grid. CPU benchmarks evaluate a BENCH_CPU_SLICE^3 voxel slice from
// the middle of it (the densest part) and are single threaded: scale by the voxel counts printed at the
// start to compare them with the GPU. The CPU density has no mdlib implementation, so the one below is
// a plain reference: per 8x8x8 block, the AOs that survive screening are evaluated once and
// rho = sum_ij D_ij phi_i phi_j is formed from them.

#define BENCH_ANGSTROM_TO_BOHR 1.8897261246257702
#define BENCH_CUTOFF   1.0e-6
#define BENCH_BLK      8
#define BENCH_BLK_VOX  (BENCH_BLK * BENCH_BLK * BENCH_BLK)
#define BENCH_CPU_SLICE 16

typedef struct {
    md_allocator_i* arena;
    md_gto_basis_t  basis;
    float*          atom_xyz;   // bohr, packed
    size_t          num_atoms;
    size_t          num_ao;
    double*         D;          // total density, [num_ao][num_ao]
    double*         mo;         // HOMO, [num_ao]
    int             mo_idx;
    md_grid_t       grid;
} bench_gto_t;

static void bench_gto_free(bench_gto_t* b) {
    if (b->arena) md_arena_allocator_destroy(b->arena);
    MEMSET(b, 0, sizeof(*b));
}

static bool bench_gto_load(bench_gto_t* b) {
    MEMSET(b, 0, sizeof(*b));
    b->arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(64));
    md_system_t sys = { .alloc = b->arena };
    md_system_state_t state = { .alloc = b->arena };
    if (!md_vlx_system_init_from_file(&sys, &state, STR_LIT(MD_BENCHMARK_DATA_DIR "/vlx/mol.h5"))) {
        MD_LOG_ERROR("Could not load benchmark vlx file");
        bench_gto_free(b);
        return false;
    }
    const md_attribute_t* c   = md_attributes_find(&sys.attributes, STR_LIT("qm/atom/coordinate"));
    const md_attribute_t* d   = md_attributes_find(&sys.attributes, STR_LIT("orbital/total/density"));
    const md_attribute_t* mo  = md_attributes_find(&sys.attributes, STR_LIT("orbital/alpha/coefficient"));
    const md_attribute_t* occ = md_attributes_find(&sys.attributes, STR_LIT("orbital/alpha/occupation"));
    if (!d) d = md_attributes_find(&sys.attributes, STR_LIT("orbital/alpha/density"));
    if (!c || !d || !mo || !occ || !md_gto_basis_extract_attributes(&b->basis, &sys.attributes, b->arena)) {
        MD_LOG_ERROR("Benchmark file lacks a basis, a geometry, a density or orbitals");
        bench_gto_free(b);
        return false;
    }

    b->num_atoms = md_attribute_value_count(&c->format);
    double* xyz = md_arena_allocator_push(b->arena, sizeof(double) * 3 * b->num_atoms);
    md_attribute_extract_f64(xyz, 3 * b->num_atoms, c, md_attribute_slice_all(), md_unit_none());
    b->atom_xyz = md_arena_allocator_push(b->arena, sizeof(float) * 3 * b->num_atoms);

    vec3_t lo = vec3_set1(FLT_MAX), hi = vec3_set1(-FLT_MAX);
    for (size_t i = 0; i < b->num_atoms; ++i) {
        const vec3_t p = vec3_set((float)(xyz[3 * i + 0] * BENCH_ANGSTROM_TO_BOHR), (float)(xyz[3 * i + 1] * BENCH_ANGSTROM_TO_BOHR), (float)(xyz[3 * i + 2] * BENCH_ANGSTROM_TO_BOHR));
        MEMCPY(b->atom_xyz + 3 * i, p.elem, sizeof(float) * 3);
        lo = vec3_min(lo, p);
        hi = vec3_max(hi, p);
    }

    const size_t n = md_gto_basis_num_ao(&b->basis);
    b->num_ao = n;
    b->D  = md_arena_allocator_push(b->arena, sizeof(double) * n * n);
    b->mo = md_arena_allocator_push(b->arena, sizeof(double) * n);
    md_attribute_extract_f64(b->D, n * n, d, md_attribute_slice_all(), md_unit_none());

    // The HOMO: the last occupied orbital
    const size_t num_mo = md_attribute_value_count(&occ->format);
    double* occupation = md_arena_allocator_push(b->arena, sizeof(double) * num_mo);
    md_attribute_extract_f64(occupation, num_mo, occ, md_attribute_slice_all(), md_unit_none());
    for (size_t i = num_mo; i-- > 0;) {
        if (occupation[i] > 0.0) { b->mo_idx = (int)i; break; }
    }
    md_attribute_extract_f64(b->mo, n, mo, md_attribute_slice_1((uint32_t)b->mo_idx), md_unit_none());

    const float h = 0.4f, margin = 6.0f;
    lo = vec3_sub1(lo, margin);
    hi = vec3_add1(hi, margin);
    b->grid = (md_grid_t){
        .orientation = mat3_ident(),
        .origin  = lo,
        .spacing = vec3_set1(h),
        .dim     = { (int)ceilf((hi.x - lo.x) / h), (int)ceilf((hi.y - lo.y) / h), (int)ceilf((hi.z - lo.z) / h) },
    };

    static bool printed = false;
    if (!printed) {
        printed = true;
        printf("  mol.h5: %zu atoms, %zu AOs, %u shells, HOMO %d, grid %dx%dx%d (%zu voxels), CPU slice %d^3 (%d voxels)\n",
               b->num_atoms, n, b->basis.num_shells, b->mo_idx + 1, b->grid.dim[0], b->grid.dim[1], b->grid.dim[2],
               md_grid_num_points(&b->grid), BENCH_CPU_SLICE, BENCH_CPU_SLICE * BENCH_CPU_SLICE * BENCH_CPU_SLICE);
    }
    return true;
}

// ---------------------------------------------------------------------------
// CPU
// ---------------------------------------------------------------------------

// The legacy CPU path evaluates at origin + index * spacing; the GPU at voxel centres. Shift the grid so
// that the CPU evaluates the same points.
static md_grid_t bench_cpu_grid(const bench_gto_t* b) {
    md_grid_t g = b->grid;
    g.origin = vec3_add(g.origin, vec3_mul1(g.spacing, 0.5f));
    return g;
}

// Index offset of the slice, a multiple of the block size
static void bench_cpu_slice_offset(int off[3], const md_grid_t* g) {
    for (int i = 0; i < 3; ++i) {
        off[i] = ((g->dim[i] - BENCH_CPU_SLICE) / 2) & ~(BENCH_BLK - 1);
    }
}

UBENCH_EX(gto, mo_cpu) {
    bench_gto_t b;
    if (!bench_gto_load(&b)) return;

    const size_t max_gtos = md_gto_pgto_count(&b.basis);
    md_gto_t* gtos     = md_arena_allocator_push(b.arena, sizeof(md_gto_t) * max_gtos);
    md_gto_t* sub_gtos = md_arena_allocator_push(b.arena, sizeof(md_gto_t) * max_gtos);
    const size_t num_gtos = md_gto_expand_with_ao_coeffs(gtos, &b.basis, b.atom_xyz, 0, b.mo, BENCH_CUTOFF);

    const md_grid_t grid = bench_cpu_grid(&b);
    float* out = md_arena_allocator_push(b.arena, sizeof(float) * md_grid_num_points(&grid));
    MEMSET(out, 0, sizeof(float) * md_grid_num_points(&grid));
    int slice_off[3];
    bench_cpu_slice_offset(slice_off, &grid);

    UBENCH_DO_BENCHMARK() {
        for (int z = 0; z < BENCH_CPU_SLICE; z += BENCH_BLK) {
            for (int y = 0; y < BENCH_CPU_SLICE; y += BENCH_BLK) {
                for (int x = 0; x < BENCH_CPU_SLICE; x += BENCH_BLK) {
                    const int off[3] = { slice_off[0] + x, slice_off[1] + y, slice_off[2] + z };
                    const int len[3] = { BENCH_BLK, BENCH_BLK, BENCH_BLK };
                    float aabb_min[3], aabb_max[3];
                    for (int i = 0; i < 3; ++i) {
                        aabb_min[i] = grid.origin.elem[i] + off[i] * grid.spacing.elem[i];
                        aabb_max[i] = grid.origin.elem[i] + (off[i] + BENCH_BLK) * grid.spacing.elem[i];
                    }
                    const size_t num_sub = md_gto_aabb_test(sub_gtos, aabb_min, aabb_max, gtos, num_gtos);
                    md_gto_grid_evaluate_sub(out, &grid, off, len, sub_gtos, num_sub, MD_GTO_EVAL_MODE_PSI);
                }
            }
        }
    }

    UBENCH_DO_NOTHING(out);
    bench_gto_free(&b);
}

// Cartesian monomial exponents and normalization of every AO of a shell, indexed [l][component]
typedef struct {
    uint8_t ijk[3];
    float   norm;
} bench_cart_t;

typedef struct {
    bench_cart_t cart[MD_GTO_MAX_ANGULAR_MOMENTUM + 1][15];
    uint32_t*    shell_ao;      // first AO of each shell
    float*       shell_r;       // radius of influence of each shell
    float*       D;             // float copy of the density matrix
    float*       phi;           // [num_ao][BENCH_BLK_VOX]
    float*       Db;            // [num_ao][num_ao], the block's sub-matrix
    uint32_t*    ao_list;
} bench_density_cpu_t;

static void bench_density_cpu_init(bench_density_cpu_t* c, const bench_gto_t* b) {
    MEMSET(c, 0, sizeof(*c));
    for (uint32_t l = 0; l <= MD_GTO_MAX_ANGULAR_MOMENTUM; ++l) {
        for (uint32_t i = 0; i < md_gto_num_cart_ao(l); ++i) {
            int x, y, z;
            md_gto_cart_ijk(&x, &y, &z, l, i);
            c->cart[l][i] = (bench_cart_t){ { (uint8_t)x, (uint8_t)y, (uint8_t)z }, (float)md_gto_cart_norm_factor(x, y, z) };
        }
    }

    const md_gto_basis_t* basis = &b->basis;
    const size_t n = b->num_ao;
    c->shell_ao = md_arena_allocator_push(b->arena, sizeof(uint32_t) * basis->num_shells);
    c->shell_r  = md_arena_allocator_push(b->arena, sizeof(float) * basis->num_shells);
    uint32_t ao = 0;
    for (uint32_t s = 0; s < basis->num_shells; ++s) {
        const md_gto_shell_t* sh = &basis->shells[s];
        c->shell_ao[s] = ao;
        ao += md_gto_num_cart_ao(sh->l);
        double r = 0.0;
        for (uint32_t p = 0; p < sh->num_primitives; ++p) {
            const uint32_t k = sh->primitive_offset + p;
            r = fmax(r, md_gto_compute_radius_of_influence((int)sh->l, 0, 0, basis->coeff[k], basis->alpha[k], BENCH_CUTOFF));
        }
        c->shell_r[s] = (float)r;
    }

    c->D = md_arena_allocator_push(b->arena, sizeof(float) * n * n);
    for (size_t i = 0; i < n * n; ++i) c->D[i] = (float)b->D[i];
    c->phi     = md_arena_allocator_push(b->arena, sizeof(float) * n * BENCH_BLK_VOX);
    c->Db      = md_arena_allocator_push(b->arena, sizeof(float) * n * n);
    c->ao_list = md_arena_allocator_push(b->arena, sizeof(uint32_t) * n);
}

static inline float bench_powi(float x, int e) {
    float r = 1.0f;
    for (; e > 0; --e) r *= x;
    return r;
}

// Evaluates the density of one 8x8x8 block of an axis aligned grid, written into the full-grid array
static void bench_density_cpu_block(float* out, const md_grid_t* grid, const int off[3], const bench_gto_t* b, const bench_density_cpu_t* c) {
    const md_gto_basis_t* basis = &b->basis;
    const size_t N = b->num_ao;

    float px[BENCH_BLK], py[BENCH_BLK], pz[BENCH_BLK];
    for (int i = 0; i < BENCH_BLK; ++i) {
        px[i] = grid->origin.x + (off[0] + i) * grid->spacing.x;
        py[i] = grid->origin.y + (off[1] + i) * grid->spacing.y;
        pz[i] = grid->origin.z + (off[2] + i) * grid->spacing.z;
    }

    // Shells whose radius of influence reaches the block, and with them the AOs of the block
    size_t n = 0;
    for (uint32_t s = 0; s < basis->num_shells; ++s) {
        const md_gto_shell_t* sh = &basis->shells[s];
        const float* ctr = b->atom_xyz + 3 * sh->atom_idx;
        const float dx = ctr[0] - fminf(fmaxf(ctr[0], px[0]), px[BENCH_BLK - 1]);
        const float dy = ctr[1] - fminf(fmaxf(ctr[1], py[0]), py[BENCH_BLK - 1]);
        const float dz = ctr[2] - fminf(fmaxf(ctr[2], pz[0]), pz[BENCH_BLK - 1]);
        const float r = c->shell_r[s];
        if (dx * dx + dy * dy + dz * dz >= r * r) continue;

        float R[BENCH_BLK_VOX], vx[BENCH_BLK_VOX], vy[BENCH_BLK_VOX], vz[BENCH_BLK_VOX];
        for (int v = 0; v < BENCH_BLK_VOX; ++v) {
            vx[v] = px[v & 7] - ctr[0];
            vy[v] = py[(v >> 3) & 7] - ctr[1];
            vz[v] = pz[v >> 6] - ctr[2];
            const float r2 = vx[v] * vx[v] + vy[v] * vy[v] + vz[v] * vz[v];
            float sum = 0.0f;
            for (uint32_t p = 0; p < sh->num_primitives; ++p) {
                sum += basis->coeff[sh->primitive_offset + p] * expf(-basis->alpha[sh->primitive_offset + p] * r2);
            }
            R[v] = sum;
        }
        for (uint32_t k = 0; k < md_gto_num_cart_ao(sh->l); ++k) {
            const bench_cart_t* ct = &c->cart[sh->l][k];
            float* row = c->phi + n * BENCH_BLK_VOX;
            for (int v = 0; v < BENCH_BLK_VOX; ++v) {
                row[v] = ct->norm * bench_powi(vx[v], ct->ijk[0]) * bench_powi(vy[v], ct->ijk[1]) * bench_powi(vz[v], ct->ijk[2]) * R[v];
            }
            c->ao_list[n++] = c->shell_ao[s] + k;
        }
    }

    float rho[BENCH_BLK_VOX] = {0};
    if (n > 0) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                c->Db[i * n + j] = c->D[(size_t)c->ao_list[i] * N + c->ao_list[j]];
            }
        }
        for (size_t i = 0; i < n; ++i) {
            float X[BENCH_BLK_VOX] = {0};
            for (size_t j = 0; j < n; ++j) {
                const float d = c->Db[i * n + j];
                const float* phi_j = c->phi + j * BENCH_BLK_VOX;
                for (int v = 0; v < BENCH_BLK_VOX; ++v) X[v] += d * phi_j[v];
            }
            const float* phi_i = c->phi + i * BENCH_BLK_VOX;
            for (int v = 0; v < BENCH_BLK_VOX; ++v) rho[v] += phi_i[v] * X[v];
        }
    }

    for (int v = 0; v < BENCH_BLK_VOX; ++v) {
        const size_t x = (size_t)off[0] + (v & 7), y = (size_t)off[1] + ((v >> 3) & 7), z = (size_t)off[2] + (v >> 6);
        out[x + (size_t)grid->dim[0] * (y + (size_t)grid->dim[1] * z)] = rho[v];
    }
}

UBENCH_EX(gto, density_cpu) {
    bench_gto_t b;
    if (!bench_gto_load(&b)) return;

    bench_density_cpu_t c;
    bench_density_cpu_init(&c, &b);

    const md_grid_t grid = bench_cpu_grid(&b);
    float* out = md_arena_allocator_push(b.arena, sizeof(float) * md_grid_num_points(&grid));
    MEMSET(out, 0, sizeof(float) * md_grid_num_points(&grid));
    int slice_off[3];
    bench_cpu_slice_offset(slice_off, &grid);

    UBENCH_DO_BENCHMARK() {
        for (int z = 0; z < BENCH_CPU_SLICE; z += BENCH_BLK) {
            for (int y = 0; y < BENCH_CPU_SLICE; y += BENCH_BLK) {
                for (int x = 0; x < BENCH_CPU_SLICE; x += BENCH_BLK) {
                    const int off[3] = { slice_off[0] + x, slice_off[1] + y, slice_off[2] + z };
                    bench_density_cpu_block(out, &grid, off, &b, &c);
                }
            }
        }
    }

    UBENCH_DO_NOTHING(out);
    bench_gto_free(&b);
}

// ---------------------------------------------------------------------------
// GPU
// ---------------------------------------------------------------------------
#if MD_ENABLE_GPU

typedef struct {
    md_gpu_device_t    device;
    md_gpu_stream_t    stream;
    md_gto_gpu_basis_t basis;
    md_gpu_addr_t      atoms;
    md_gpu_addr_t      coeff;   // sized for the density, the larger of the two packings
    md_gpu_texture_t   tex;
} bench_gpu_t;

static void bench_gpu_free(bench_gpu_t* g) {
    if (g->tex)   md_gpu_texture_destroy(g->tex);
    if (g->coeff) md_gpu_free(g->stream, g->coeff);
    if (g->atoms) md_gpu_free(g->stream, g->atoms);
    if (g->basis) md_gto_gpu_basis_destroy(g->stream, g->basis);
    if (g->stream) md_gpu_stream_sync(g->stream);
    if (g->device) {
        md_gto_gpu_shutdown();
        md_gpu_device_destroy(g->device);
    }
    MEMSET(g, 0, sizeof(*g));
}

// False, with nothing left allocated, when there is no usable device
static bool bench_gpu_init(bench_gpu_t* g, const bench_gto_t* b) {
    MEMSET(g, 0, sizeof(*g));
    g->device = md_gpu_device_create(NULL);
    if (!g->device) return false;
    g->stream = md_gpu_stream_default(g->device, MD_GPU_STREAM_COMPUTE);
    md_gto_gpu_initialize(g->device);

    const md_gto_gpu_basis_desc_t basis_desc = { .basis = &b->basis, .cutoff = BENCH_CUTOFF };
    g->basis = md_gto_gpu_basis_create(g->stream, &basis_desc);
    const size_t num_cgtos = md_gto_gpu_basis_num_cgtos(g->basis);

    const size_t atoms_size = md_gto_gpu_atom_buffer_size(b->num_atoms);
    g->atoms = md_gpu_malloc(g->stream, MD_GPU_MEM_DEVICE, atoms_size).gpu;
    g->coeff = md_gpu_malloc(g->stream, MD_GPU_MEM_DEVICE, md_gto_gpu_coeff_size_density(num_cgtos)).gpu;
    g->tex = md_gpu_texture_create(g->stream, &(md_gpu_texture_desc_t){
        .type = MD_GPU_TEX_3D, .format = MD_GPU_FORMAT_R32_FLOAT, .usage = MD_GPU_TEX_STORAGE,
        .width = (uint32_t)b->grid.dim[0], .height = (uint32_t)b->grid.dim[1], .depth_or_layers = (uint32_t)b->grid.dim[2],
    });
    if (!g->basis || num_cgtos != b->num_ao || !g->atoms || !g->coeff || !g->tex) {
        MD_LOG_ERROR("Failed to set up the GPU benchmark");
        bench_gpu_free(g);
        return false;
    }

    float* dst = md_gpu_upload_begin(g->stream, g->atoms, atoms_size);
    md_gto_gpu_atom_pack(dst, b->atom_xyz, 0, b->num_atoms);
    md_gpu_upload_end(g->stream);
    return true;
}

static void bench_gpu_upload_mo(bench_gpu_t* g, const bench_gto_t* b) {
    const double* mo[1] = { b->mo };
    float* dst = md_gpu_upload_begin(g->stream, g->coeff, md_gto_gpu_coeff_size_mo(1, b->num_ao));
    md_gto_gpu_coeff_pack_mo(dst, mo, NULL, 1, b->num_ao);
    md_gpu_upload_end(g->stream);
}

static void bench_gpu_upload_density(bench_gpu_t* g, const bench_gto_t* b) {
    float* dst = md_gpu_upload_begin(g->stream, g->coeff, md_gto_gpu_coeff_size_density(b->num_ao));
    md_gto_gpu_coeff_pack_density(dst, b->D, b->num_ao);
    md_gpu_upload_end(g->stream);
}

#define BENCH_GPU_MO(NAME, ALGO) \
UBENCH_EX(gto, NAME) { \
    bench_gto_t b; \
    if (!bench_gto_load(&b)) return; \
    bench_gpu_t g; \
    if (!bench_gpu_init(&g, &b)) { bench_gto_free(&b); return; } \
    bench_gpu_upload_mo(&g, &b); \
    const md_gto_gpu_orbital_desc_t desc = { \
        .basis = g.basis, .atom_xyz = g.atoms, .coeff = g.coeff, .out_tex = g.tex, .grid = &b.grid, \
        .sample_offset = {0.5f, 0.5f, 0.5f}, .num_orbitals = 1, .eval_mode = MD_GTO_EVAL_MODE_PSI, \
        .op = MD_GTO_OP_SET, .algorithm = ALGO, \
    }; \
    /* Warm up: kernel creation, first-use compilation */ \
    md_gto_gpu_orbital_launch(g.stream, &desc); \
    md_gpu_stream_sync(g.stream); \
    UBENCH_DO_BENCHMARK() { \
        md_gto_gpu_orbital_launch(g.stream, &desc); \
        md_gpu_stream_sync(g.stream); \
    } \
    bench_gpu_free(&g); \
    bench_gto_free(&b); \
}

#define BENCH_GPU_DENSITY(NAME, ALGO) \
UBENCH_EX(gto, NAME) { \
    bench_gto_t b; \
    if (!bench_gto_load(&b)) return; \
    bench_gpu_t g; \
    if (!bench_gpu_init(&g, &b)) { bench_gto_free(&b); return; } \
    bench_gpu_upload_density(&g, &b); \
    const md_gto_gpu_density_desc_t desc = { \
        .basis = g.basis, .atom_xyz = g.atoms, .coeff = g.coeff, .out_tex = g.tex, .grid = &b.grid, \
        .sample_offset = {0.5f, 0.5f, 0.5f}, .op = MD_GTO_OP_SET, .algorithm = ALGO, \
    }; \
    md_gto_gpu_density_launch(g.stream, &desc); \
    md_gpu_stream_sync(g.stream); \
    UBENCH_DO_BENCHMARK() { \
        md_gto_gpu_density_launch(g.stream, &desc); \
        md_gpu_stream_sync(g.stream); \
    } \
    bench_gpu_free(&g); \
    bench_gto_free(&b); \
}

BENCH_GPU_MO(mo_gpu_reference, MD_GTO_GPU_ORBITAL_ALGO_REFERENCE)
BENCH_GPU_MO(mo_gpu_shell,     MD_GTO_GPU_ORBITAL_ALGO_SHELL)
BENCH_GPU_MO(mo_gpu_gemm,      MD_GTO_GPU_ORBITAL_ALGO_GEMM)

BENCH_GPU_DENSITY(density_gpu_reference, MD_GTO_GPU_DENSITY_ALGO_REFERENCE)
BENCH_GPU_DENSITY(density_gpu_gemm,      MD_GTO_GPU_DENSITY_ALGO_GEMM)

#endif

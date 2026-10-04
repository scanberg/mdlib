// bench_gto_gpu.c  (md_bench_gto_gpu)
//
// Correctness and performance comparison of the GPU GTO kernels, electron densities
// (md_gto_gpu_density_launch) and molecular orbitals (md_gto_gpu_orbital_launch), meant
// to be run on every GPU we care about. The log starts with a description of the system
// (OS, CPU, memory, compiler, build, GPU and driver) so results stay attributable.
//
//   md_bench_gto_gpu [--quick] [--iters N] [--seconds S] [--case NAME]... [--dim N]
//                    [--algos LIST] [--scratch-mb MB] [--data DIR]
//                    [--device SEL] [--prefer high-performance|low-power] [--list-devices]
//
// GPU selection: --device takes an index from --list-devices or a case-insensitive part
// of the adapter name ("intel", "1060"); --prefer ranks the adapters otherwise (and among
// several matches). Without either, MD_GPU_DEVICE is honoured, then the discrete GPU.
//
// For every test case and grid size, each algorithm is
//   * run once (warm-up, includes any runtime shader compilation) and checked: against
//     the reference kernel over the whole grid, and against a double-precision CPU
//     evaluation without screening at 2048 sampled voxels,
//   * then timed (wall clock, launch to stream sync) until ~1 s or the iteration cap
//     has been spent; the median is reported.
//
// Density algorithms (--algos takes a comma separated list; default: all listed here):
//   reference                 original per-AO kernel
//   tiled                     shell-based register-tiled kernel
//   gemm-v1                   first two-pass GEMM kernel (64 voxels x 64 AOs per group)
//   gemm                      two-pass GEMM, configuration chosen from the GPU vendor
//   gemm-<V>x<T>[-s<N>]       two-pass GEMM with V voxels per group (64/128/256), T-wide
//                             AO tiles (32/64), and 32-wide tiles for blocks with <= N AOs
//   gemm-otf[-<V>x<T>...]     as gemm, gathering D inside the GEMM pass
//   gemm-sg, gemm-sg-otf      simdgroup_matrix GEMM (Metal only; elsewhere = gemm-v1)
//
// Orbital algorithms, run for 1 orbital (psi), 1 orbital (psi^2) and 32 orbitals (psi^2):
//   mo-reference              original per-AO kernel
//   mo                        automatic (shell kernel up to 8 orbitals, GEMM above)
//   mo-shell[-v<N>][-exact]   shell kernel, N voxels per thread (1/2/4; default per
//                             vendor), -exact: no coefficient-aware screening
//   mo-gemm[-<T>]             two-pass Phi + GEMM, T orbitals per tile (32/64)
//
// Test cases (Cartesian AOs, cutoff 1e-6 as used by VIAMD). Geometries and the def2-SVP
// basis come from benchmark/data/density, so no test_data checkout and no HDF5 are needed:
//   mol    26-atom organic molecule (C, N, H), def2-SVP. With HDF5 and test_data present,
//          the SCF density and orbitals of test_data/vlx/mol.h5 are used; otherwise
//          synthetic ones.
//   c60    C60, def2-SVP, synthetic symmetric D and orbitals
//   c60f   as c60, plus a d and an f polarisation shell per atom (higher l)
//   c240   240 atoms of C720, def2-SVP, synthetic D and orbitals

#include <md_gto.h>
#include <md_system.h>
#include <md_attributes.h>
#include <core/md_gpu.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_log.h>
#include <core/md_str.h>
#include <core/md_vec_math.h>
#include <core/md_grid.h>
#ifdef MD_HDF5
#include <md_vlx.h>
#endif

#if defined(_WIN32)
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#  if defined(_M_X64) || defined(_M_IX86)
#    include <intrin.h>
#  endif
#  ifdef _MSC_VER
#    pragma comment(lib, "advapi32.lib")
#  endif
#elif defined(__APPLE__)
#  include <sys/types.h>
#  include <sys/sysctl.h>
#  include <sys/utsname.h>
#  include <unistd.h>
#else
#  include <sys/utsname.h>
#  include <unistd.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

#define ANG_TO_BOHR 1.8897261246257702
#define CUTOFF 1.0e-6
#define NUM_SAMPLES 2048
#define MAX_ALGOS 64
#define BENCH_NUM_MOS 32     // orbitals in the many-orbital test

#ifndef MD_DENSITY_DATA_DIR
#define MD_DENSITY_DATA_DIR "data/density"
#endif

// ---------------------------------------------------------------------------
// Algorithms
// ---------------------------------------------------------------------------

typedef enum { ALGO_DENSITY, ALGO_ORBITAL } algo_kind_t;

typedef struct {
    char name[48];
    algo_kind_t kind;
    // density
    md_gto_gpu_density_algo_t algo;
    uint32_t gp, gm;
    int32_t  sn;
    // orbitals
    md_gto_gpu_orbital_algo_t oalgo;
    uint32_t vpt;
    bool     exact;
    uint32_t otile;
} algo_spec_t;

static const char* DEFAULT_ALGOS =
    "reference,tiled,gemm-v1,gemm,"
    "gemm-64x32,gemm-64x64,gemm-64x64-s32,"
    "gemm-128x32,gemm-128x64,gemm-128x64-s32,"
    "gemm-256x32,gemm-256x64,gemm-256x64-s32,"
    "mo-reference,mo,mo-shell-v1,mo-shell-v2,mo-shell-v4,mo-shell-v2-exact,mo-gemm-32,mo-gemm-64";

static bool parse_algo(algo_spec_t* out, const char* name) {
    memset(out, 0, sizeof(*out));
    snprintf(out->name, sizeof(out->name), "%s", name);
    if (!strncmp(name, "mo", 2)) {
        out->kind = ALGO_ORBITAL;
        const char* rest = name + 2;
        if (!strcmp(rest, ""))           { out->oalgo = MD_GTO_GPU_ORBITAL_ALGO_DEFAULT;   return true; }
        if (!strcmp(rest, "-reference")) { out->oalgo = MD_GTO_GPU_ORBITAL_ALGO_REFERENCE; return true; }
        if (!strncmp(rest, "-shell", 6)) {
            out->oalgo = MD_GTO_GPU_ORBITAL_ALGO_SHELL;
            rest += 6;
            unsigned v = 0;
            int used = 0;
            if (sscanf(rest, "-v%u%n", &v, &used) == 1) {
                if (v != 1 && v != 2 && v != 4) return false;
                out->vpt = v;
                rest += used;
            }
            if (!strcmp(rest, "-exact")) { out->exact = true; rest += 6; }
            return *rest == 0;
        }
        if (!strncmp(rest, "-gemm", 5)) {
            out->oalgo = MD_GTO_GPU_ORBITAL_ALGO_GEMM;
            rest += 5;
            if (*rest == 0) return true;
            unsigned t = 0;
            if (sscanf(rest, "-%u", &t) != 1 || (t != 32 && t != 64)) return false;
            out->otile = t;
            return true;
        }
        return false;
    }
    out->kind = ALGO_DENSITY;
    if (!strcmp(name, "reference"))   { out->algo = MD_GTO_GPU_DENSITY_ALGO_REFERENCE;   return true; }
    if (!strcmp(name, "tiled"))       { out->algo = MD_GTO_GPU_DENSITY_ALGO_TILED;       return true; }
    if (!strcmp(name, "gemm-v1"))     { out->algo = MD_GTO_GPU_DENSITY_ALGO_GEMM_V1;     return true; }
    if (!strcmp(name, "gemm-sg"))     { out->algo = MD_GTO_GPU_DENSITY_ALGO_GEMM_SG;     return true; }
    if (!strcmp(name, "gemm-sg-otf")) { out->algo = MD_GTO_GPU_DENSITY_ALGO_GEMM_SG_OTF; return true; }

    const char* rest = NULL;
    if (!strncmp(name, "gemm-otf", 8)) { out->algo = MD_GTO_GPU_DENSITY_ALGO_GEMM_OTF; rest = name + 8; }
    else if (!strncmp(name, "gemm", 4)) { out->algo = MD_GTO_GPU_DENSITY_ALGO_GEMM;    rest = name + 4; }
    else return false;
    if (*rest == 0) return true;   // automatic configuration
    unsigned gp = 0, gm = 0;
    int sn = 0, used = 0;
    if (sscanf(rest, "-%ux%u%n", &gp, &gm, &used) != 2) return false;
    rest += used;
    if (*rest) {
        if (sscanf(rest, "-s%d", &sn) != 1) return false;
    } else {
        sn = -1;   // explicit configuration without a small-block split
    }
    if ((gp != 64 && gp != 128 && gp != 256) || (gm != 32 && gm != 64)) return false;
    out->gp = gp; out->gm = gm; out->sn = sn;
    return true;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static uint64_t rng_state = 0x9E3779B97F4A7C15ull;
static double rng_uniform(void) {   // [0,1)
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (double)(rng_state >> 11) * (1.0 / 9007199254740992.0);
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}

static bool streq_nocase(const char* a, const char* b) {
    for (; *a && *b; ++a, ++b) {
        if (toupper((unsigned char)*a) != toupper((unsigned char)*b)) return false;
    }
    return *a == *b;
}

typedef struct {
    float* xyz;           // bohr
    char (*element)[4];
    size_t count;
} geometry_t;

static bool read_xyz(geometry_t* g, const char* path) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "could not open %s\n", path); return false; }
    char line[512];
    size_t n = 0;
    if (!fgets(line, sizeof(line), f) || sscanf(line, "%zu", &n) != 1 || !fgets(line, sizeof(line), f)) { fclose(f); return false; }
    g->xyz = (float*)malloc(sizeof(float) * 3 * n);
    g->element = malloc(sizeof(*g->element) * n);
    size_t i = 0;
    while (i < n && fgets(line, sizeof(line), f)) {
        char el[16];
        double x, y, z;
        if (sscanf(line, "%15s %lf %lf %lf", el, &x, &y, &z) == 4) {
            snprintf(g->element[i], sizeof(g->element[i]), "%s", el);
            g->xyz[3 * i + 0] = (float)(x * ANG_TO_BOHR);
            g->xyz[3 * i + 1] = (float)(y * ANG_TO_BOHR);
            g->xyz[3 * i + 2] = (float)(z * ANG_TO_BOHR);
            i++;
        }
    }
    fclose(f);
    g->count = i;
    return i == n;
}

// ---------------------------------------------------------------------------
// Basis sets in VeloxChem format, normalised to the md_gto_basis_t convention:
// coeff carries the radial normalisation of the x^l component of the shell.
// ---------------------------------------------------------------------------

#define MAX_EL_SHELLS 32
#define MAX_EL_PRIMS  128

typedef struct {
    char     element[4];
    uint32_t num_shells;
    uint32_t l[MAX_EL_SHELLS];
    uint32_t num_prims[MAX_EL_SHELLS];
    uint32_t prim_offset[MAX_EL_SHELLS];
    double   alpha[MAX_EL_PRIMS];
    double   coeff[MAX_EL_PRIMS];
} element_basis_t;

static double dfact_odd(int l) {   // (2l-1)!!
    double r = 1.0;
    for (int k = 2 * l - 1; k > 1; k -= 2) r *= k;
    return r;
}

static double prim_norm(double alpha, int l) {
    return pow(2.0 * alpha / 3.14159265358979323846, 0.75) * pow(4.0 * alpha, 0.5 * l) / sqrt(dfact_odd(l));
}

// Same normalisation as the VeloxChem reader (md_vlx.c, normalize_basis_set), so the
// basis matches what VIAMD builds from a VeloxChem file with this basis set.
static void normalise_shell(double* coeff, const double* alpha, int n, int l) {
    static const double F[] = { 0.0, 2.0, 1.15470053837925152902, 1.03279555898864450271, 0.19518001458970663587 };
    static const double OV[] = { 1.0, 0.5, 3.0, 7.5, 420.0 };
    const double pi = 3.14159265358979323846;
    if (n == 1) coeff[0] = 1.0;
    for (int i = 0; i < n; ++i) {
        coeff[i] *= pow(alpha[i] * 2.0 / pi, 0.75) * pow(F[l] * alpha[i], 0.5 * l);
    }
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double fab = 1.0 / (alpha[i] + alpha[j]);
            s += coeff[i] * coeff[j] * pow(pi * fab, 1.5) * OV[l] * pow(fab, l);
        }
    }
    const double inv = 1.0 / sqrt(s);
    for (int i = 0; i < n; ++i) coeff[i] *= inv;
}

static int read_basis_file(element_basis_t* out, int cap, const char* path) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "could not open %s\n", path); return 0; }
    char line[512];
    int count = 0;
    element_basis_t* cur = NULL;
    while (fgets(line, sizeof(line), f)) {
        char tag[32], el[8];
        if (sscanf(line, "%31s", tag) != 1) continue;
        if (!strcmp(tag, "@ATOMBASIS") && sscanf(line, "%31s %7s", tag, el) == 2 && count < cap) {
            cur = &out[count++];
            memset(cur, 0, sizeof(*cur));
            snprintf(cur->element, sizeof(cur->element), "%s", el);
            continue;
        }
        if (!strcmp(tag, "@END")) { cur = NULL; continue; }
        if (!cur) continue;
        char lch;
        int np = 0, nc = 0;
        if (sscanf(line, " %c %d %d", &lch, &np, &nc) == 3 && isalpha((unsigned char)lch)) {
            static const char LS[] = "SPDFG";
            const char* lp = strchr(LS, toupper((unsigned char)lch));
            if (!lp || np <= 0 || cur->num_shells >= MAX_EL_SHELLS) continue;
            const uint32_t sh = cur->num_shells++;
            const uint32_t off = sh ? cur->prim_offset[sh - 1] + cur->num_prims[sh - 1] : 0;
            cur->l[sh] = (uint32_t)(lp - LS);
            cur->num_prims[sh] = (uint32_t)np;
            cur->prim_offset[sh] = off;
            for (int p = 0; p < np && off + p < MAX_EL_PRIMS; ++p) {
                if (!fgets(line, sizeof(line), f)) break;
                sscanf(line, "%lf %lf", &cur->alpha[off + p], &cur->coeff[off + p]);
            }
            normalise_shell(&cur->coeff[off], &cur->alpha[off], np, (int)cur->l[sh]);
        }
    }
    fclose(f);
    return count;
}

// Basis for a geometry from per-element bases, plus optional extra single-primitive
// shells (l, alpha) on every atom.
static bool build_basis(md_gto_basis_t* out, const geometry_t* g, const element_basis_t* eb, int num_eb,
                        const int* extra_l, const double* extra_alpha, int num_extra) {
    uint32_t ns = 0, np = 0;
    const element_basis_t** per_atom = malloc(sizeof(*per_atom) * g->count);
    for (size_t a = 0; a < g->count; ++a) {
        per_atom[a] = NULL;
        for (int e = 0; e < num_eb; ++e) {
            if (streq_nocase(eb[e].element, g->element[a])) { per_atom[a] = &eb[e]; break; }
        }
        if (!per_atom[a]) { fprintf(stderr, "no basis for element %s\n", g->element[a]); free(per_atom); return false; }
        ns += per_atom[a]->num_shells + (uint32_t)num_extra;
        for (uint32_t s = 0; s < per_atom[a]->num_shells; ++s) np += per_atom[a]->num_prims[s];
        np += (uint32_t)num_extra;
    }
    out->num_shells = ns;
    out->num_primitives = np;
    out->shells = (md_gto_shell_t*)malloc(sizeof(md_gto_shell_t) * ns);
    out->alpha  = (float*)malloc(sizeof(float) * np);
    out->coeff  = (float*)malloc(sizeof(float) * np);
    uint32_t si = 0, pi = 0;
    for (size_t a = 0; a < g->count; ++a) {
        const element_basis_t* e = per_atom[a];
        for (uint32_t s = 0; s < e->num_shells; ++s) {
            out->shells[si++] = (md_gto_shell_t){ .atom_idx = (uint32_t)a, .primitive_offset = pi, .num_primitives = e->num_prims[s], .l = e->l[s] };
            for (uint32_t p = 0; p < e->num_prims[s]; ++p, ++pi) {
                out->alpha[pi] = (float)e->alpha[e->prim_offset[s] + p];
                out->coeff[pi] = (float)e->coeff[e->prim_offset[s] + p];
            }
        }
        for (int x = 0; x < num_extra; ++x) {
            out->shells[si++] = (md_gto_shell_t){ .atom_idx = (uint32_t)a, .primitive_offset = pi, .num_primitives = 1, .l = (uint32_t)extra_l[x] };
            out->alpha[pi] = (float)extra_alpha[x];
            out->coeff[pi] = (float)prim_norm(extra_alpha[x], extra_l[x]);
            pi++;
        }
    }
    free(per_atom);
    return true;
}

typedef struct {
    char            name[32];
    md_gto_basis_t  basis;       // heap allocated
    float*          atom_xyz;    // bohr, 3 per atom
    size_t          num_atoms;
    double*         D;           // full num_ao x num_ao, symmetric
    size_t          num_ao;
    const char*     d_kind;
    double*         C;           // BENCH_NUM_MOS orbitals x num_ao; row 0 is the single-orbital test
    const char*     c_kind;
} bench_case_t;

static double* ao_centers(const md_gto_basis_t* basis, const float* xyz, size_t num_ao) {
    double* c = (double*)malloc(sizeof(double) * 3 * num_ao);
    size_t k = 0;
    for (uint32_t s = 0; s < basis->num_shells; ++s) {
        uint32_t n = md_gto_num_cart_ao(basis->shells[s].l);
        for (uint32_t i = 0; i < n; ++i, ++k) {
            const float* p = xyz + 3 * basis->shells[s].atom_idx;
            c[3 * k + 0] = p[0]; c[3 * k + 1] = p[1]; c[3 * k + 2] = p[2];
        }
    }
    return c;
}

// Synthetic orbitals: random coefficients localised around an atom (decay length 8 bohr),
// a different centre per orbital, so that coefficient-aware screening sees realistic
// variation between orbitals.
static double* synthetic_C(const md_gto_basis_t* basis, const float* xyz, size_t num_atoms, size_t num_ao) {
    double* c = ao_centers(basis, xyz, num_ao);
    double* C = (double*)malloc(sizeof(double) * BENCH_NUM_MOS * num_ao);
    for (size_t m = 0; m < BENCH_NUM_MOS; ++m) {
        const float* ctr = xyz + 3 * ((m * 7919) % num_atoms);
        for (size_t i = 0; i < num_ao; ++i) {
            double dx = c[3*i+0] - ctr[0], dy = c[3*i+1] - ctr[1], dz = c[3*i+2] - ctr[2];
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            C[m * num_ao + i] = (2.0 * rng_uniform() - 1.0) * 0.3 * exp(-r / 8.0);
        }
    }
    free(c);
    return C;
}

// Synthetic symmetric density-like matrix: random entries damped with the distance
// between the AO centres, positive diagonal.
static double* synthetic_D(const md_gto_basis_t* basis, const float* xyz, size_t num_ao) {
    double* c = ao_centers(basis, xyz, num_ao);
    double* D = (double*)malloc(sizeof(double) * num_ao * num_ao);
    for (size_t i = 0; i < num_ao; ++i) {
        D[i * num_ao + i] = 0.2 + 0.8 * rng_uniform();
        for (size_t j = i + 1; j < num_ao; ++j) {
            double dx = c[3*i+0] - c[3*j+0], dy = c[3*i+1] - c[3*j+1], dz = c[3*i+2] - c[3*j+2];
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            double v = (2.0 * rng_uniform() - 1.0) * 0.5 * exp(-0.25 * r);
            D[i * num_ao + j] = v;
            D[j * num_ao + i] = v;
        }
    }
    free(c);
    return D;
}

static bool case_from_files(bench_case_t* bc, const char* name, const char* data_dir, const char* xyz_file,
                            const element_basis_t* eb, int num_eb, const int* el, const double* ea, int ne) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/%s", data_dir, xyz_file);
    geometry_t g = {0};
    if (!read_xyz(&g, path)) return false;
    bool ok = build_basis(&bc->basis, &g, eb, num_eb, el, ea, ne);
    free(g.element);
    if (!ok) { free(g.xyz); return false; }
    bc->atom_xyz  = g.xyz;
    bc->num_atoms = g.count;
    bc->num_ao    = md_gto_basis_num_ao(&bc->basis);
    bc->D         = synthetic_D(&bc->basis, bc->atom_xyz, bc->num_ao);
    bc->C         = synthetic_C(&bc->basis, bc->atom_xyz, bc->num_atoms, bc->num_ao);
    snprintf(bc->name, sizeof(bc->name), "%s", name);
    bc->d_kind = "synthetic";
    bc->c_kind = "synthetic";
    return true;
}

#ifdef MD_HDF5
// The mol case with the SCF density of test_data/vlx/mol.h5, when available.
static bool case_mol_scf(bench_case_t* bc) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(64));
    md_system_t sys = { .alloc = arena };
    md_system_state_t state = { .alloc = arena };
    bool ok = false;
    FILE* probe = fopen(MD_BENCHMARK_DATA_DIR "/vlx/mol.h5", "rb");
    if (!probe) goto done;
    fclose(probe);
    if (!md_vlx_system_init_from_file(&sys, &state, STR_LIT(MD_BENCHMARK_DATA_DIR "/vlx/mol.h5"))) goto done;

    md_gto_basis_t basis = {0};
    if (!md_gto_basis_extract_attributes(&basis, &sys.attributes, arena)) goto done;
    const md_attribute_t* coord_attr = md_attributes_find(&sys.attributes, STR_LIT("qm/atom/coordinate"));
    const md_attribute_t* c_attr = md_attributes_find(&sys.attributes, STR_LIT("orbital/alpha/coefficient"));
    const md_attribute_t* o_attr = md_attributes_find(&sys.attributes, STR_LIT("orbital/alpha/occupation"));
    if (!coord_attr || !c_attr || !o_attr) goto done;

    const size_t nv = md_attribute_element_count(&coord_attr->format);
    double* xyz = (double*)malloc(sizeof(double) * nv);
    md_attribute_extract_f64(xyz, nv, coord_attr, md_attribute_slice_all(), md_unit_none());
    bc->num_atoms = nv / 3;
    bc->atom_xyz = (float*)malloc(sizeof(float) * nv);
    for (size_t i = 0; i < nv; ++i) bc->atom_xyz[i] = (float)(xyz[i] * ANG_TO_BOHR);
    free(xyz);

    bc->basis = basis;
    bc->basis.shells = (md_gto_shell_t*)malloc(sizeof(md_gto_shell_t) * basis.num_shells);
    bc->basis.alpha  = (float*)malloc(sizeof(float) * basis.num_primitives);
    bc->basis.coeff  = (float*)malloc(sizeof(float) * basis.num_primitives);
    memcpy(bc->basis.shells, basis.shells, sizeof(md_gto_shell_t) * basis.num_shells);
    memcpy(bc->basis.alpha,  basis.alpha,  sizeof(float) * basis.num_primitives);
    memcpy(bc->basis.coeff,  basis.coeff,  sizeof(float) * basis.num_primitives);
    bc->num_ao = md_gto_basis_num_ao(&bc->basis);

    // D = sum_i n_i c_i c_i^T (closed shell: n = 2).
    const size_t num_mo = md_attribute_element_count(&o_attr->format);
    double* occ = (double*)malloc(sizeof(double) * num_mo);
    md_attribute_extract_f64(occ, num_mo, o_attr, md_attribute_slice_all(), md_unit_none());
    const size_t n = bc->num_ao;
    double* coeff = (double*)malloc(sizeof(double) * n);
    bc->D = (double*)calloc(n * n, sizeof(double));
    double max_occ = 0.0;
    for (size_t i = 0; i < num_mo; ++i) max_occ = fmax(max_occ, occ[i]);
    const double scale = (max_occ <= 1.0) ? 2.0 : 1.0;
    for (size_t i = 0; i < num_mo; ++i) {
        if (occ[i] == 0.0) continue;
        if (md_attribute_extract_f64(coeff, n, c_attr, md_attribute_slice_1((uint32_t)i), md_unit_none()) != n) continue;
        const double w = scale * occ[i];
        for (size_t a = 0; a < n; ++a) {
            const double wa = w * coeff[a];
            for (size_t b = 0; b < n; ++b) bc->D[a * n + b] += wa * coeff[b];
        }
    }
    // Orbitals: HOMO first, then downwards (the single-orbital test uses the HOMO).
    size_t homo = 0;
    for (size_t i = 0; i < num_mo; ++i) if (occ[i] != 0.0) homo = i;
    bc->C = (double*)calloc(BENCH_NUM_MOS * n, sizeof(double));
    for (size_t m = 0; m < BENCH_NUM_MOS && m <= homo; ++m) {
        md_attribute_extract_f64(bc->C + m * n, n, c_attr, md_attribute_slice_1((uint32_t)(homo - m)), md_unit_none());
    }
    free(coeff);
    free(occ);
    snprintf(bc->name, sizeof(bc->name), "mol");
    bc->d_kind = "SCF (occupied MOs of mol.h5)";
    bc->c_kind = "SCF (HOMO and below, mol.h5)";
    ok = true;
done:
    md_arena_allocator_destroy(arena);
    return ok;
}
#endif

// ---------------------------------------------------------------------------
// CPU reference
// ---------------------------------------------------------------------------

// All Cartesian AO values at point p (bohr), double precision, no screening.
static void cpu_eval_aos(double* out, const md_gto_basis_t* basis, const float* xyz, const double p[3]) {
    size_t k = 0;
    for (uint32_t s = 0; s < basis->num_shells; ++s) {
        const md_gto_shell_t* sh = &basis->shells[s];
        const float* c = xyz + 3 * sh->atom_idx;
        const double dx = p[0] - c[0], dy = p[1] - c[1], dz = p[2] - c[2];
        const double r2 = dx*dx + dy*dy + dz*dz;
        double R = 0.0;
        for (uint32_t i = 0; i < sh->num_primitives; ++i) {
            R += (double)basis->coeff[sh->primitive_offset + i] * exp(-(double)basis->alpha[sh->primitive_offset + i] * r2);
        }
        const uint32_t ncart = md_gto_num_cart_ao(sh->l);
        for (uint32_t ci = 0; ci < ncart; ++ci, ++k) {
            int i, j, l;
            md_gto_cart_ijk(&i, &j, &l, sh->l, ci);
            out[k] = md_gto_cart_norm_factor(i, j, l) * pow(dx, i) * pow(dy, j) * pow(dz, l) * R;
        }
    }
}

static double cpu_density(const bench_case_t* bc, const double p[3], double* phi, uint32_t* idx) {
    cpu_eval_aos(phi, &bc->basis, bc->atom_xyz, p);
    // Products below 1e-14 relative to O(1) densities are irrelevant; keep the rest.
    size_t m = 0;
    for (size_t k = 0; k < bc->num_ao; ++k) {
        if (fabs(phi[k]) > 1.0e-12) idx[m++] = (uint32_t)k;
    }
    const size_t n = bc->num_ao;
    double rho = 0.0;
    for (size_t a = 0; a < m; ++a) {
        const double* Drow = bc->D + (size_t)idx[a] * n;
        double t = 0.0;
        for (size_t b = 0; b < m; ++b) t += Drow[idx[b]] * phi[idx[b]];
        rho += phi[idx[a]] * t;
    }
    return rho;
}

// sum_m f(psi_m) over the first `count` orbitals of bc->C, f = psi or psi^2.
static double cpu_orbitals(const bench_case_t* bc, const double p[3], double* phi, uint32_t* idx, size_t count, bool squared) {
    cpu_eval_aos(phi, &bc->basis, bc->atom_xyz, p);
    size_t m = 0;
    for (size_t k = 0; k < bc->num_ao; ++k) {
        if (fabs(phi[k]) > 1.0e-12) idx[m++] = (uint32_t)k;
    }
    double res = 0.0;
    for (size_t o = 0; o < count; ++o) {
        const double* c = bc->C + o * bc->num_ao;
        double psi = 0.0;
        for (size_t a = 0; a < m; ++a) psi += c[idx[a]] * phi[idx[a]];
        res += squared ? psi * psi : psi;
    }
    return res;
}

// World position of voxel (i,j,k) for an axis aligned grid sampled at voxel centres.
static void voxel_pos(double out[3], const md_grid_t* g, int i, int j, int k) {
    out[0] = g->origin.x + (i + 0.5) * g->spacing.x;
    out[1] = g->origin.y + (j + 0.5) * g->spacing.y;
    out[2] = g->origin.z + (k + 0.5) * g->spacing.z;
}

// Average / max number of screened AOs per 8^3 block, using the same shell radii
// (max over components and primitives) as the GPU kernels.
static void block_stats(double* out_avg, uint32_t* out_max, double* out_nonempty, const bench_case_t* bc, const md_grid_t* g) {
    const md_gto_basis_t* basis = &bc->basis;
    double* r = (double*)malloc(sizeof(double) * basis->num_shells);
    for (uint32_t s = 0; s < basis->num_shells; ++s) {
        const md_gto_shell_t* sh = &basis->shells[s];
        double max_r = 0.0;
        for (uint32_t ci = 0; ci < md_gto_num_cart_ao(sh->l); ++ci) {
            int i, j, l;
            md_gto_cart_ijk(&i, &j, &l, sh->l, ci);
            const double nrm = md_gto_cart_norm_factor(i, j, l);
            for (uint32_t p = 0; p < sh->num_primitives; ++p) {
                const float coeff = (float)(basis->coeff[sh->primitive_offset + p] * nrm);
                max_r = fmax(max_r, md_gto_compute_radius_of_influence(i, j, l, coeff, basis->alpha[sh->primitive_offset + p], CUTOFF));
            }
        }
        r[s] = max_r;
    }
    const int nb[3] = { (g->dim[0] + 7) / 8, (g->dim[1] + 7) / 8, (g->dim[2] + 7) / 8 };
    uint64_t sum = 0, nonempty = 0;
    uint32_t mx = 0;
    for (int bz = 0; bz < nb[2]; ++bz)
    for (int by = 0; by < nb[1]; ++by)
    for (int bx = 0; bx < nb[0]; ++bx) {
        const double mn[3] = { g->origin.x + bx * 8 * g->spacing.x, g->origin.y + by * 8 * g->spacing.y, g->origin.z + bz * 8 * g->spacing.z };
        const double mxp[3] = { mn[0] + 8 * g->spacing.x, mn[1] + 8 * g->spacing.y, mn[2] + 8 * g->spacing.z };
        uint32_t cnt = 0;
        for (uint32_t s = 0; s < basis->num_shells; ++s) {
            if (r[s] <= 0.0) continue;
            const float* c = bc->atom_xyz + 3 * basis->shells[s].atom_idx;
            double d2 = 0.0;
            for (int a = 0; a < 3; ++a) {
                double v = c[a] < mn[a] ? mn[a] - c[a] : (c[a] > mxp[a] ? c[a] - mxp[a] : 0.0);
                d2 += v * v;
            }
            if (d2 < r[s] * r[s]) cnt += md_gto_num_cart_ao(basis->shells[s].l);
        }
        sum += cnt;
        nonempty += (cnt > 0);
        if (cnt > mx) mx = cnt;
    }
    const double total = (double)nb[0] * nb[1] * nb[2];
    *out_avg = nonempty ? (double)sum / (double)nonempty : 0.0;
    *out_max = mx;
    *out_nonempty = nonempty / total;
    free(r);
}

// ---------------------------------------------------------------------------
// GPU driver
// ---------------------------------------------------------------------------

static size_t g_scratch_bytes = 0;

typedef struct {
    md_gpu_device_t  dev;
    md_gpu_stream_t  stream;
} gpu_ctx_t;

typedef struct {
    double ms_median;
    double ms_min;
    int    iters;
} algo_result_t;

typedef struct {
    md_gpu_addr_t      coeff;      // density matrix (density) or orbital coefficients (orbitals)
    size_t             num_mos;    // orbitals
    md_gto_eval_mode_t mode;       // orbitals
} run_input_t;

static void launch(gpu_ctx_t* ctx, md_gto_gpu_basis_t gb, md_gpu_addr_t atoms, md_gpu_texture_t tex, const md_grid_t* grid,
                   const algo_spec_t* spec, const run_input_t* in) {
    if (spec->kind == ALGO_DENSITY) {
        md_gto_gpu_density_desc_t desc = {
            .basis = gb, .atom_xyz = atoms, .coeff = in->coeff, .out_tex = tex, .grid = grid,
            .sample_offset = {0.5f, 0.5f, 0.5f}, .op = MD_GTO_OP_SET, .algorithm = spec->algo,
            .scratch_bytes = g_scratch_bytes,
            .gemm_group_voxels = spec->gp, .gemm_tile = spec->gm, .gemm_small_block = spec->sn,
        };
        md_gto_gpu_density_launch(ctx->stream, &desc);
    } else {
        md_gto_gpu_orbital_desc_t desc = {
            .basis = gb, .atom_xyz = atoms, .coeff = in->coeff, .out_tex = tex, .grid = grid,
            .sample_offset = {0.5f, 0.5f, 0.5f}, .num_orbitals = in->num_mos, .eval_mode = in->mode,
            .op = MD_GTO_OP_SET, .algorithm = spec->oalgo, .voxels_per_thread = spec->vpt,
            .exact_screening = spec->exact, .gemm_tile = spec->otile, .scratch_bytes = g_scratch_bytes,
        };
        md_gto_gpu_orbital_launch(ctx->stream, &desc);
    }
}

static void run_algo(algo_result_t* res, float* out_grid, gpu_ctx_t* ctx, md_gto_gpu_basis_t gb, md_gpu_addr_t atoms,
                     md_gpu_texture_t tex, md_gpu_mem_t readback, const md_grid_t* grid, const algo_spec_t* spec,
                     const run_input_t* in, int max_iters, double min_seconds) {
    const size_t num_vox = (size_t)grid->dim[0] * grid->dim[1] * grid->dim[2];

    // Clear, so that a kernel that silently fails to run cannot pass on stale data.
    float* zero = (float*)calloc(num_vox, sizeof(float));
    md_gpu_upload_texture(ctx->stream, tex, NULL, zero, sizeof(float) * num_vox);
    free(zero);

    // Warm-up (shader compilation, allocator warm-up); also the run that is checked.
    launch(ctx, gb, atoms, tex, grid, spec, in);
    md_gpu_copy_from_texture(ctx->stream, readback.gpu, tex, NULL);
    md_gpu_stream_sync(ctx->stream);
    memcpy(out_grid, readback.cpu, sizeof(float) * num_vox);

    double times[64];
    int n = 0;
    double total = 0.0;
    while (n < max_iters && n < 64 && (n < 2 || total < min_seconds)) {
        md_tick_t t0 = md_tick_now();
        launch(ctx, gb, atoms, tex, grid, spec, in);
        md_gpu_stream_sync(ctx->stream);
        md_tick_t t1 = md_tick_now();
        times[n] = md_tick_to_milliseconds(t1 - t0);
        total += times[n] / 1000.0;
        n++;
    }
    qsort(times, n, sizeof(double), cmp_double);
    res->ms_median = times[n / 2];
    res->ms_min = times[0];
    res->iters = n;
}

typedef struct {
    gpu_ctx_t*         ctx;
    const bench_case_t* bc;
    md_gto_gpu_basis_t gb;
    md_gpu_addr_t      atoms;
    md_gpu_texture_t   tex;
    md_gpu_mem_t       readback;
    const md_grid_t*   grid;
    float              mn[3];
    float              h;
    float*             ref_grid;
    float*             cur_grid;
    double*            phi;
    uint32_t*          idx;
    int                max_iters;
    double             min_seconds;
} table_ctx_t;

// One table: every algorithm of `kind`, the reference one first.
static void run_table(table_ctx_t* T, const char* title, algo_kind_t kind, const run_input_t* in,
                      const algo_spec_t* algos, int num_algos) {
    const bench_case_t* bc = T->bc;
    const int* gd = T->grid->dim;
    const size_t num_vox = (size_t)gd[0] * gd[1] * gd[2];

    int num_kind = 0;
    for (int a = 0; a < num_algos; ++a) num_kind += algos[a].kind == kind;
    if (num_kind == 0) return;

    // CPU reference at sampled voxels: half uniform, half near atoms.
    int    sample_idx[3 * NUM_SAMPLES];
    double* sample_ref = (double*)malloc(sizeof(double) * NUM_SAMPLES);
    double ref_max = 0.0;
    for (int s = 0; s < NUM_SAMPLES; ++s) {
        int i, j, k;
        if (s & 1) {
            const float* c = bc->atom_xyz + 3 * (size_t)(rng_uniform() * bc->num_atoms);
            i = (int)((c[0] - T->mn[0]) / T->h + (rng_uniform() - 0.5) * 6.0);
            j = (int)((c[1] - T->mn[1]) / T->h + (rng_uniform() - 0.5) * 6.0);
            k = (int)((c[2] - T->mn[2]) / T->h + (rng_uniform() - 0.5) * 6.0);
            i = i < 0 ? 0 : (i >= gd[0] ? gd[0] - 1 : i);
            j = j < 0 ? 0 : (j >= gd[1] ? gd[1] - 1 : j);
            k = k < 0 ? 0 : (k >= gd[2] ? gd[2] - 1 : k);
        } else {
            i = (int)(rng_uniform() * gd[0]);
            j = (int)(rng_uniform() * gd[1]);
            k = (int)(rng_uniform() * gd[2]);
        }
        sample_idx[3 * s + 0] = i; sample_idx[3 * s + 1] = j; sample_idx[3 * s + 2] = k;
        double p[3];
        voxel_pos(p, T->grid, i, j, k);
        sample_ref[s] = (kind == ALGO_DENSITY)
            ? cpu_density(bc, p, T->phi, T->idx)
            : cpu_orbitals(bc, p, T->phi, T->idx, in->num_mos, in->mode == MD_GTO_EVAL_MODE_PSI_SQUARED);
        ref_max = fmax(ref_max, fabs(sample_ref[s]));
    }

    printf("\n  [%s]\n", title);
    printf("  %-20s %10s %10s %8s %6s %12s %12s\n", "algorithm", "median ms", "min ms", "speedup", "iters", "err vs CPU", "diff vs ref");
    double ref_ms = 0.0, ref_gmax = 0.0;
    bool have_ref = false;
    double best_ms = DBL_MAX;
    const char* best_name = "-";
    const md_gto_gpu_density_algo_t ref_d = MD_GTO_GPU_DENSITY_ALGO_REFERENCE;
    const md_gto_gpu_orbital_algo_t ref_o = MD_GTO_GPU_ORBITAL_ALGO_REFERENCE;
    // Reference first, then the rest in list order.
    for (int pass = 0; pass < 2; ++pass) {
        for (int a = 0; a < num_algos; ++a) {
            if (algos[a].kind != kind) continue;
            const bool is_ref = (kind == ALGO_DENSITY) ? algos[a].algo == ref_d : algos[a].oalgo == ref_o;
            if (is_ref != (pass == 0)) continue;

            algo_result_t r = {0};
            float* g = is_ref ? T->ref_grid : T->cur_grid;
            run_algo(&r, g, T->ctx, T->gb, T->atoms, T->tex, T->readback, T->grid, &algos[a], in, T->max_iters, T->min_seconds);

            double err = 0.0;
            for (int s = 0; s < NUM_SAMPLES; ++s) {
                const int i = sample_idx[3*s+0], j = sample_idx[3*s+1], k = sample_idx[3*s+2];
                err = fmax(err, fabs(g[((size_t)k * gd[1] + j) * gd[0] + i] - sample_ref[s]));
            }
            if (is_ref) {
                ref_ms = r.ms_median;
                have_ref = true;
                for (size_t v = 0; v < num_vox; ++v) ref_gmax = fmax(ref_gmax, fabs((double)g[v]));
            }
            char diffs[32] = "-";
            if (have_ref) {
                double diff = 0.0;
                if (!is_ref) for (size_t v = 0; v < num_vox; ++v) diff = fmax(diff, fabs((double)g[v] - (double)T->ref_grid[v]));
                snprintf(diffs, sizeof(diffs), "%.2e", diff / fmax(ref_gmax, 1e-30));
            }
            char speedup[16] = "-";
            if (ref_ms > 0.0) snprintf(speedup, sizeof(speedup), "%.2fx", ref_ms / r.ms_median);
            const double rel_err = err / fmax(ref_max, 1e-30);
            printf("  %-20s %10.2f %10.2f %8s %6d %12.2e %12s%s\n", algos[a].name, r.ms_median, r.ms_min,
                   speedup, r.iters, rel_err, diffs, rel_err > 1.0e-3 ? "  <-- WRONG RESULT" : "");
            if (rel_err <= 1.0e-3 && r.ms_median < best_ms) { best_ms = r.ms_median; best_name = algos[a].name; }
        }
    }
    printf("  (errors relative to max |value| = %.4g; CPU reference without screening)\n", ref_max);
    if (best_ms < DBL_MAX) {
        printf("  fastest: %s, %.2f ms, %.2fx vs reference\n", best_name, best_ms, ref_ms > 0.0 ? ref_ms / best_ms : 0.0);
    }
    free(sample_ref);
}

static void run_case(gpu_ctx_t* ctx, const bench_case_t* bc, const int* dims, int num_dims, int max_iters, double min_seconds,
                     const algo_spec_t* algos, int num_algos) {
    md_gto_gpu_basis_t gb = md_gto_gpu_basis_create(ctx->stream, &(md_gto_gpu_basis_desc_t){ .basis = &bc->basis, .cutoff = CUTOFF });
    if (!gb) { printf("  basis upload failed\n"); return; }

    const size_t n = bc->num_ao;
    md_gpu_addr_t atoms  = md_gpu_malloc(ctx->stream, MD_GPU_MEM_DEVICE, md_gto_gpu_atom_buffer_size(bc->num_atoms)).gpu;
    md_gpu_addr_t dmat   = md_gpu_malloc(ctx->stream, MD_GPU_MEM_DEVICE, md_gto_gpu_coeff_size_density(n)).gpu;
    md_gpu_addr_t mocoef = md_gpu_malloc(ctx->stream, MD_GPU_MEM_DEVICE, md_gto_gpu_coeff_size_mo(BENCH_NUM_MOS, n)).gpu;
    {
        float* p = (float*)md_gpu_upload_begin(ctx->stream, atoms, md_gto_gpu_atom_buffer_size(bc->num_atoms));
        md_gto_gpu_atom_pack(p, bc->atom_xyz, 0, bc->num_atoms);
        md_gpu_upload_end(ctx->stream);
        p = (float*)md_gpu_upload_begin(ctx->stream, dmat, md_gto_gpu_coeff_size_density(n));
        md_gto_gpu_coeff_pack_density(p, bc->D, n);
        md_gpu_upload_end(ctx->stream);
        const double* rows[BENCH_NUM_MOS];
        for (int m = 0; m < BENCH_NUM_MOS; ++m) rows[m] = bc->C + (size_t)m * n;
        p = (float*)md_gpu_upload_begin(ctx->stream, mocoef, md_gto_gpu_coeff_size_mo(BENCH_NUM_MOS, n));
        md_gto_gpu_coeff_pack_mo(p, rows, NULL, BENCH_NUM_MOS, n);
        md_gpu_upload_end(ctx->stream);
    }

    // Bounding box + margin.
    float mn[3] = { FLT_MAX, FLT_MAX, FLT_MAX }, mx[3] = { -FLT_MAX, -FLT_MAX, -FLT_MAX };
    for (size_t i = 0; i < bc->num_atoms; ++i) {
        for (int a = 0; a < 3; ++a) {
            mn[a] = fminf(mn[a], bc->atom_xyz[3 * i + a]);
            mx[a] = fmaxf(mx[a], bc->atom_xyz[3 * i + a]);
        }
    }
    const float margin = 6.0f; // bohr
    for (int a = 0; a < 3; ++a) { mn[a] -= margin; mx[a] += margin; }

    printf("\n=== case %s: %zu atoms, %u shells, %zu Cartesian AOs, D %s, orbitals %s ===\n",
           bc->name, bc->num_atoms, bc->basis.num_shells, n, bc->d_kind, bc->c_kind);

    double* phi = (double*)malloc(sizeof(double) * n);
    uint32_t* idx = (uint32_t*)malloc(sizeof(uint32_t) * n);

    for (int di = 0; di < num_dims; ++di) {
        const int dim = dims[di];
        // Cubic voxels, longest side gets 'dim' voxels.
        float ext = fmaxf(mx[0] - mn[0], fmaxf(mx[1] - mn[1], mx[2] - mn[2]));
        float h = ext / (float)dim;
        int gd[3];
        for (int a = 0; a < 3; ++a) gd[a] = (int)ceilf((mx[a] - mn[a]) / h);
        md_grid_t grid = {
            .orientation = mat3_ident(),
            .origin = vec3_set(mn[0], mn[1], mn[2]),
            .spacing = vec3_set(h, h, h),
            .dim = { gd[0], gd[1], gd[2] },
        };
        const size_t num_vox = (size_t)gd[0] * gd[1] * gd[2];

        double avg_ao; uint32_t max_ao; double frac_nonempty;
        block_stats(&avg_ao, &max_ao, &frac_nonempty, bc, &grid);
        printf("\n-- grid %d x %d x %d (%.3f bohr), screened AOs per non-empty 8^3 block: avg %.0f, max %u, %.0f%% of blocks non-empty\n",
               gd[0], gd[1], gd[2], h, avg_ao, max_ao, 100.0 * frac_nonempty);

        table_ctx_t T = {
            .ctx = ctx, .bc = bc, .gb = gb, .atoms = atoms, .grid = &grid,
            .mn = { mn[0], mn[1], mn[2] }, .h = h, .phi = phi, .idx = idx,
            .max_iters = max_iters, .min_seconds = min_seconds,
        };
        T.tex = md_gpu_texture_create(ctx->stream, &(md_gpu_texture_desc_t){
            .type = MD_GPU_TEX_3D, .format = MD_GPU_FORMAT_R32_FLOAT, .usage = MD_GPU_TEX_STORAGE,
            .width = (uint32_t)gd[0], .height = (uint32_t)gd[1], .depth_or_layers = (uint32_t)gd[2],
        });
        T.readback = md_gpu_malloc(ctx->stream, MD_GPU_MEM_HOST_READ, sizeof(float) * num_vox);
        T.ref_grid = (float*)malloc(sizeof(float) * num_vox);
        T.cur_grid = (float*)malloc(sizeof(float) * num_vox);
        if (!T.tex || !T.readback.cpu || !T.ref_grid || !T.cur_grid) { printf("  allocation failed\n"); return; }

        run_input_t din = { .coeff = dmat };
        run_table(&T, "density", ALGO_DENSITY, &din, algos, num_algos);

        run_input_t o1  = { .coeff = mocoef, .num_mos = 1, .mode = MD_GTO_EVAL_MODE_PSI };
        run_table(&T, "1 orbital, psi", ALGO_ORBITAL, &o1, algos, num_algos);
        run_input_t o1s = { .coeff = mocoef, .num_mos = 1, .mode = MD_GTO_EVAL_MODE_PSI_SQUARED };
        run_table(&T, "1 orbital, psi^2", ALGO_ORBITAL, &o1s, algos, num_algos);
        run_input_t o32 = { .coeff = mocoef, .num_mos = BENCH_NUM_MOS, .mode = MD_GTO_EVAL_MODE_PSI_SQUARED };
        run_table(&T, "32 orbitals, sum of psi^2", ALGO_ORBITAL, &o32, algos, num_algos);

        free(T.ref_grid);
        free(T.cur_grid);
        md_gpu_free(ctx->stream, T.readback.gpu);
        md_gpu_texture_destroy(T.tex);
    }
    free(phi);
    free(idx);
    md_gpu_free(ctx->stream, mocoef);
    md_gpu_free(ctx->stream, dmat);
    md_gpu_free(ctx->stream, atoms);
    md_gto_gpu_basis_destroy(ctx->stream, gb);
    md_gpu_stream_sync(ctx->stream);
}

// ---------------------------------------------------------------------------
// System description for the log
// ---------------------------------------------------------------------------

#ifndef MD_BENCH_GIT_REV
#define MD_BENCH_GIT_REV "unknown"
#endif
#ifndef MD_BENCH_BUILD_TYPE
#define MD_BENCH_BUILD_TYPE ""
#endif
#ifdef NDEBUG
#define MD_BENCH_NDEBUG "NDEBUG"
#else
#define MD_BENCH_NDEBUG "no NDEBUG"
#endif

static void trim(char* s) {
    size_t n = strlen(s);
    while (n && isspace((unsigned char)s[n - 1])) s[--n] = 0;
    size_t b = 0;
    while (s[b] && isspace((unsigned char)s[b])) ++b;
    if (b) memmove(s, s + b, n - b + 1);
}

static void compiler_string(char* out, size_t cap) {
#if defined(__clang__)
    snprintf(out, cap, "clang %s", __clang_version__);
#elif defined(_MSC_VER)
    snprintf(out, cap, "MSVC %d (_MSC_FULL_VER %d)", _MSC_VER, _MSC_FULL_VER);
#elif defined(__GNUC__)
    snprintf(out, cap, "gcc %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#else
    snprintf(out, cap, "unknown");
#endif
}

#if defined(_WIN32)
static bool reg_string(const char* key, const char* value, char* out, DWORD cap) {
    out[0] = 0;
    return RegGetValueA(HKEY_LOCAL_MACHINE, key, value, RRF_RT_REG_SZ, NULL, out, &cap) == ERROR_SUCCESS;
}
static DWORD reg_dword(const char* key, const char* value) {
    DWORD v = 0, cap = sizeof(v);
    if (RegGetValueA(HKEY_LOCAL_MACHINE, key, value, RRF_RT_REG_DWORD, NULL, &v, &cap) != ERROR_SUCCESS) return 0;
    return v;
}
#elif defined(__APPLE__)
static bool sysctl_string(const char* name, char* out, size_t cap) {
    out[0] = 0;
    size_t len = cap;
    if (sysctlbyname(name, out, &len, NULL, 0) != 0) { out[0] = 0; return false; }
    out[cap - 1] = 0;
    return true;
}
static int64_t sysctl_int(const char* name) {
    int64_t v = 0;
    size_t len = sizeof(v);
    if (sysctlbyname(name, &v, &len, NULL, 0) != 0) return -1;
    if (len == sizeof(int32_t)) { int32_t w; memcpy(&w, &v, sizeof(w)); return w; }
    return v;
}
#else
// First "key<sep>value" line of a text file whose key matches; value trimmed.
static bool file_field(const char* path, const char* key, char sep, char* out, size_t cap) {
    out[0] = 0;
    FILE* f = fopen(path, "r");
    if (!f) return false;
    char line[512];
    const size_t klen = strlen(key);
    bool found = false;
    while (fgets(line, sizeof(line), f)) {
        if (strncmp(line, key, klen)) continue;
        const char* p = line + klen;
        while (*p == ' ' || *p == '\t') ++p;
        if (*p != sep) continue;
        snprintf(out, cap, "%s", p + 1);
        trim(out);
        size_t n = strlen(out);   // strip quotes (os-release)
        if (n >= 2 && out[0] == '"' && out[n - 1] == '"') { memmove(out, out + 1, n - 2); out[n - 2] = 0; }
        found = true;
        break;
    }
    fclose(f);
    return found;
}
#endif

static const char* device_type_name(md_gpu_device_type_t t) {
    switch (t) {
    case MD_GPU_DEVICE_TYPE_DISCRETE:   return "discrete";
    case MD_GPU_DEVICE_TYPE_INTEGRATED: return "integrated";
    case MD_GPU_DEVICE_TYPE_VIRTUAL:    return "virtual";
    case MD_GPU_DEVICE_TYPE_CPU:        return "cpu";
    default:                            return "other";
    }
}

static void print_system_info(int argc, char** argv, const md_gpu_device_info_t* gpu, const char* gpu_vendor) {
    char os[512] = "unknown", cpu[256] = "unknown", cores[128] = "unknown", host[256] = "unknown", comp[256];
    double ram_gib = 0.0;

#if defined(_WIN32)
    {
        const char* cv = "SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion";
        char product[128], display[64], build[32];
        reg_string(cv, "ProductName", product, sizeof(product));
        reg_string(cv, "DisplayVersion", display, sizeof(display));
        reg_string(cv, "CurrentBuildNumber", build, sizeof(build));
        const DWORD ubr = reg_dword(cv, "UBR");
        // ProductName still says "Windows 10" on Windows 11; the build number tells (>= 22000).
        const bool win11 = atoi(build) >= 22000;
        if (win11 && !strncmp(product, "Windows 10", 10)) { product[8] = '1'; product[9] = '1'; }
        snprintf(os, sizeof(os), "%s %s (build %s.%lu)", product, display, build, (unsigned long)ubr);
        reg_string("HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0", "ProcessorNameString", cpu, sizeof(cpu));
        trim(cpu);
        SYSTEM_INFO si;
        GetNativeSystemInfo(&si);
        DWORD len = 0;
        int physical = 0;
        GetLogicalProcessorInformation(NULL, &len);
        SYSTEM_LOGICAL_PROCESSOR_INFORMATION* lpi = (SYSTEM_LOGICAL_PROCESSOR_INFORMATION*)malloc(len ? len : 1);
        if (lpi && GetLogicalProcessorInformation(lpi, &len)) {
            for (DWORD i = 0; i < len / sizeof(*lpi); ++i) physical += lpi[i].Relationship == RelationProcessorCore;
        }
        free(lpi);
        snprintf(cores, sizeof(cores), "%d cores, %lu logical processors", physical, (unsigned long)si.dwNumberOfProcessors);
        MEMORYSTATUSEX ms = { .dwLength = sizeof(ms) };
        if (GlobalMemoryStatusEx(&ms)) ram_gib = (double)ms.ullTotalPhys / (1024.0 * 1024.0 * 1024.0);
        DWORD hl = sizeof(host);
        if (!GetComputerNameA(host, &hl)) snprintf(host, sizeof(host), "unknown");
    }
#elif defined(__APPLE__)
    {
        char ver[64], build[64], model[128];
        struct utsname u;
        sysctl_string("kern.osproductversion", ver, sizeof(ver));
        sysctl_string("kern.osversion", build, sizeof(build));
        sysctl_string("hw.model", model, sizeof(model));
        if (uname(&u) == 0) snprintf(os, sizeof(os), "macOS %s (%s), %s %s %s, model %s", ver, build, u.sysname, u.release, u.machine, model);
        sysctl_string("machdep.cpu.brand_string", cpu, sizeof(cpu));
        const int64_t phys = sysctl_int("hw.physicalcpu"), logi = sysctl_int("hw.logicalcpu");
        const int64_t p0 = sysctl_int("hw.perflevel0.physicalcpu"), p1 = sysctl_int("hw.perflevel1.physicalcpu");
        if (p0 > 0 && p1 > 0) snprintf(cores, sizeof(cores), "%lld cores (%lld performance + %lld efficiency), %lld logical",
                                       (long long)phys, (long long)p0, (long long)p1, (long long)logi);
        else snprintf(cores, sizeof(cores), "%lld cores, %lld logical", (long long)phys, (long long)logi);
        const int64_t mem = sysctl_int("hw.memsize");
        if (mem > 0) ram_gib = (double)mem / (1024.0 * 1024.0 * 1024.0);
        gethostname(host, sizeof(host));
    }
#else
    {
        char pretty[256];
        struct utsname u;
        if (!file_field("/etc/os-release", "PRETTY_NAME", '=', pretty, sizeof(pretty))) snprintf(pretty, sizeof(pretty), "Linux");
        if (uname(&u) == 0) snprintf(os, sizeof(os), "%s, kernel %s %s", pretty, u.release, u.machine);
        if (!file_field("/proc/cpuinfo", "model name", ':', cpu, sizeof(cpu)) &&
            !file_field("/proc/cpuinfo", "Model", ':', cpu, sizeof(cpu))) {
            snprintf(cpu, sizeof(cpu), "unknown");
        }
        snprintf(cores, sizeof(cores), "%ld logical processors online", sysconf(_SC_NPROCESSORS_ONLN));
        char mt[64];
        if (file_field("/proc/meminfo", "MemTotal", ':', mt, sizeof(mt))) ram_gib = atof(mt) / (1024.0 * 1024.0);
        gethostname(host, sizeof(host));
    }
#endif
    host[sizeof(host) - 1] = 0;
    compiler_string(comp, sizeof(comp));

    char date[64] = "unknown";
    time_t now = time(NULL);
    struct tm* tmv = gmtime(&now);
    if (tmv) strftime(date, sizeof(date), "%Y-%m-%d %H:%M:%S UTC", tmv);

    printf("# ---------------------------------------------------------------------------\n");
    printf("# md_bench_gto_gpu\n");
    printf("# date:      %s\n", date);
    printf("# host:      %s\n", host);
    printf("# os:        %s\n", os);
    printf("# cpu:       %s, %s\n", cpu, cores);
    printf("# memory:    %.1f GiB\n", ram_gib);
    printf("# compiler:  %s, %zu-bit\n", comp, sizeof(void*) * 8);
    const char* build_type = MD_BENCH_BUILD_TYPE;
    printf("# build:     %s (%s), mdlib %s\n", build_type[0] ? build_type : "no build type", MD_BENCH_NDEBUG, MD_BENCH_GIT_REV);
    printf("# gpu:       %s (vendor 0x%04X %s, device 0x%04X, %s)\n", gpu->name, gpu->vendor_id, gpu_vendor, gpu->device_id,
           device_type_name(gpu->type));
    printf("# gpu info:  subgroup %u, max %u threads per group\n", gpu->preferred_group_multiple, gpu->max_threads_per_group);
    printf("# driver:    %s\n", gpu->driver[0] ? gpu->driver : "unknown");
    md_gpu_adapter_info_t adapters[16];
    const uint32_t num_adapters = md_gpu_enumerate_adapters(adapters, 16);
    for (uint32_t i = 0; i < num_adapters && i < 16; ++i) {
        const md_gpu_adapter_info_t* a = &adapters[i];
        printf("# adapter %u: %s (0x%04X:0x%04X, %s)%s%s%s\n", i, a->name, a->vendor_id, a->device_id, device_type_name(a->type),
               a->usable ? "" : " unusable, lacks ", a->usable ? "" : a->missing,
               (gpu->name[0] && i == gpu->adapter_index && !strcmp(a->name, gpu->name)) ? "  <-- used" : "");
    }
    const char* sel = getenv("MD_GPU_DEVICE");
    if (sel && sel[0]) printf("# MD_GPU_DEVICE=%s\n", sel);
    printf("# command:  ");
    for (int i = 0; i < argc; ++i) printf(" %s", argv[i]);
    printf("\n# ---------------------------------------------------------------------------\n");
}

static const char* vendor_name(uint32_t id) {
    switch (id) {
    case 0x10DE: return "NVIDIA";
    case 0x1002: return "AMD";
    case 0x8086: return "Intel";
    case 0x106B: return "Apple";
    case 0x10005: return "Mesa (software)";
    default: return "unknown";
    }
}

int main(int argc, char** argv) {
#ifdef _WIN32
    setvbuf(stdout, NULL, _IONBF, 0);   // keep stdout and the library's log lines in order
#else
    setvbuf(stdout, NULL, _IOLBF, 0);
#endif
    bool quick = false;
    int max_iters = 10;
    double min_seconds = 1.0;
    const char* only_cases[8] = {0};
    int num_only = 0;
    int dim_override = 0;
    const char* algo_list = DEFAULT_ALGOS;
    const char* data_dir = MD_DENSITY_DATA_DIR;
    const char* device_sel = NULL;
    md_gpu_device_preference_t preference = MD_GPU_DEVICE_PREFER_DEFAULT;
    bool list_devices = false;

    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--quick")) { quick = true; }
        else if (!strcmp(argv[i], "--iters") && i + 1 < argc) { max_iters = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--seconds") && i + 1 < argc) { min_seconds = atof(argv[++i]); }
        else if (!strcmp(argv[i], "--dim") && i + 1 < argc) { dim_override = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--scratch-mb") && i + 1 < argc) { g_scratch_bytes = (size_t)atoi(argv[++i]) << 20; }
        else if (!strcmp(argv[i], "--case") && i + 1 < argc && num_only < 8) { only_cases[num_only++] = argv[++i]; }
        else if (!strcmp(argv[i], "--algos") && i + 1 < argc) { algo_list = argv[++i]; }
        else if (!strcmp(argv[i], "--data") && i + 1 < argc) { data_dir = argv[++i]; }
        else if (!strcmp(argv[i], "--device") && i + 1 < argc) { device_sel = argv[++i]; }
        else if (!strcmp(argv[i], "--list-devices")) { list_devices = true; }
        else if (!strcmp(argv[i], "--prefer") && i + 1 < argc) {
            const char* v = argv[++i];
            if (!strcmp(v, "high-performance") || !strcmp(v, "discrete"))  preference = MD_GPU_DEVICE_PREFER_HIGH_PERFORMANCE;
            else if (!strcmp(v, "low-power") || !strcmp(v, "integrated")) preference = MD_GPU_DEVICE_PREFER_LOW_POWER;
            else { fprintf(stderr, "--prefer takes high-performance or low-power\n"); return 1; }
        }
        else { fprintf(stderr, "unknown argument %s\n", argv[i]); return 1; }
    }
    if (quick) { max_iters = 2; min_seconds = 0.0; }

    // Algorithms. The reference kernel of each kind that appears in the list is always
    // run too (first, as the baseline of its tables).
    algo_spec_t algos[MAX_ALGOS];
    int num_algos = 2;
    bool have_kind[2] = { false, false };
    {
        char buf[2048];
        snprintf(buf, sizeof(buf), "%s", algo_list);
        for (char* tok = strtok(buf, ", "); tok && num_algos < MAX_ALGOS; tok = strtok(NULL, ", ")) {
            algo_spec_t a;
            if (!parse_algo(&a, tok)) { fprintf(stderr, "unknown algorithm '%s'\n", tok); return 1; }
            have_kind[a.kind] = true;
            if (!strcmp(tok, "reference") || !strcmp(tok, "mo-reference")) continue;
            algos[num_algos++] = a;
        }
    }
    {
        // Drop the reference slot of a kind that is not in the list.
        int k = 0;
        algo_spec_t refs[2];
        parse_algo(&refs[0], "reference");
        parse_algo(&refs[1], "mo-reference");
        if (have_kind[ALGO_DENSITY]) algos[k++] = refs[0];
        if (have_kind[ALGO_ORBITAL]) algos[k++] = refs[1];
        if (k < 2) memmove(algos + k, algos + 2, sizeof(algo_spec_t) * (size_t)(num_algos - 2));
        num_algos -= 2 - k;
    }

    if (list_devices) {
        md_gpu_adapter_info_t adapters[16];
        const uint32_t n = md_gpu_enumerate_adapters(adapters, 16);
        if (n == 0) printf("No GPU adapters: %s\n", md_gpu_last_error() ? md_gpu_last_error() : "none found");
        for (uint32_t i = 0; i < n && i < 16; ++i) {
            const md_gpu_adapter_info_t* a = &adapters[i];
            printf("%u: %s (vendor 0x%04X, device 0x%04X, %s)%s%s\n   driver: %s\n", i, a->name, a->vendor_id, a->device_id,
                   device_type_name(a->type), a->usable ? "" : " UNUSABLE, lacks ", a->usable ? "" : a->missing, a->driver);
        }
        return 0;
    }

    gpu_ctx_t ctx = {0};
    ctx.dev = md_gpu_device_create(&(md_gpu_device_desc_t){ .adapter = device_sel, .preference = preference, .label = "md_bench_gto_gpu" });
    md_gpu_device_info_t info = {0};
    if (!ctx.dev) {
        char err[2560];   // copied first: the system description enumerates adapters, which resets it
        snprintf(err, sizeof(err), "%s", md_gpu_last_error() ? md_gpu_last_error() : "unknown error");
        snprintf(info.name, sizeof(info.name), "none");
        print_system_info(argc, argv, &info, "-");
        printf("No GPU device: %s\n", err);
        return 1;
    }
    md_gpu_device_info(ctx.dev, &info);
    print_system_info(argc, argv, &info, vendor_name(info.vendor_id));
    printf("# algorithms:");
    for (int a = 0; a < num_algos; ++a) printf(" %s", algos[a].name);
    printf("\n");
    ctx.stream = md_gpu_stream_default(ctx.dev, MD_GPU_STREAM_COMPUTE);
    md_gto_gpu_initialize(ctx.dev);

    element_basis_t eb[8];
    char path[1024];
    snprintf(path, sizeof(path), "%s/def2-svp.txt", data_dir);
    const int num_eb = read_basis_file(eb, 8, path);
    if (num_eb == 0) {
        fprintf(stderr, "Could not read the basis set from %s (use --data DIR)\n", path);
        return 1;
    }

    typedef struct { const char* name; int dims[3]; int num_dims; } case_cfg_t;
    case_cfg_t cfgs[] = {
        {"mol",  {128, 256}, 2},
        {"c60",  {128, 256}, 2},
        {"c60f", {128},      1},
        {"c240", {128},      1},
    };
    for (size_t c = 0; c < sizeof(cfgs) / sizeof(cfgs[0]); ++c) {
        if (num_only) {
            bool found = false;
            for (int k = 0; k < num_only; ++k) found |= !strcmp(only_cases[k], cfgs[c].name);
            if (!found) continue;
        }
        bench_case_t bc = {0};
        bool ok = false;
        if (!strcmp(cfgs[c].name, "mol")) {
#ifdef MD_HDF5
            ok = case_mol_scf(&bc);
#endif
            if (!ok) ok = case_from_files(&bc, "mol", data_dir, "mol.xyz", eb, num_eb, NULL, NULL, 0);
        } else if (!strcmp(cfgs[c].name, "c60")) {
            ok = case_from_files(&bc, "c60", data_dir, "c60.xyz", eb, num_eb, NULL, NULL, 0);
        } else if (!strcmp(cfgs[c].name, "c60f")) {
            const int el[2] = {2, 3};
            const double ea[2] = {0.55, 0.80};
            ok = case_from_files(&bc, "c60f", data_dir, "c60.xyz", eb, num_eb, el, ea, 2);
        } else if (!strcmp(cfgs[c].name, "c240")) {
            ok = case_from_files(&bc, "c240", data_dir, "c240.xyz", eb, num_eb, NULL, NULL, 0);
        }
        if (!ok) { printf("case %s: setup failed\n", cfgs[c].name); continue; }
        int dims[3];
        int nd = cfgs[c].num_dims;
        for (int k = 0; k < nd; ++k) dims[k] = cfgs[c].dims[k];
        if (dim_override) { dims[0] = dim_override; nd = 1; }
        if (quick) { dims[0] = MIN(dims[0], 64); nd = 1; }
        run_case(&ctx, &bc, dims, nd, max_iters, min_seconds, algos, num_algos);

        free(bc.basis.shells); free(bc.basis.alpha); free(bc.basis.coeff);
        free(bc.atom_xyz); free(bc.D); free(bc.C);
    }

    md_gto_gpu_shutdown();
    md_gpu_device_destroy(ctx.dev);
    return 0;
}

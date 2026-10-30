// md_gto_iso_gpu.inl - the md_gpu path of md_gto_iso (included by md_gto_iso.c under MD_ENABLE_GPU).
//
// Host side of src/shaders/iso/gto_iso_gpu.slang: builds a field set (basis, coefficient matrix C,
// terms per field), then issues the pipeline as a job (md_gto_iso_gpu_job_*): octree level by level
// (bricks in Morton order, their corner slots), corner values (for densities whose maxima are at the
// atoms, from seeds and a flood), Surface Nets, relaxation on the lattice's surface, Newton
// projection, per-field regrouping. The host waits for the device three times (the brick and slot
// counts, the vertex count, the final counts with the mesh), and once per flood round.
// md_gto_iso_gpu_evaluate_launch runs the same vertex evaluators on the caller's points (at the end).

#include <core/md_intrinsics.h>
#include <core/md_gpu.h>
#include <gto_iso_gpu_shaders.inl>

// BENCHMARK (temporary, branch gto-iso-gpu-opt): which round of changes this is, in every profile report
#define ISOG_REVISION "r37 (the vertex evaluators at arbitrary points: md_gto_iso_gpu_evaluate; r36: the D form in tiles of the AO pairs that meet, built on the device)"

typedef struct isog_shell_t {        // mirrors IsoShell
    md_gpu_float4 c;
    md_gpu_uint4  info;
    md_gpu_uint4  info2;
} isog_shell_t;

typedef struct isog_field_t {        // mirrors IsoField
    float    iso;
    float    slack;
    uint32_t term_off;
    uint32_t term_count;
} isog_field_t;

typedef struct isog_term_t {         // mirrors IsoTerm
    float    w;
    uint32_t a, b, _pad;
} isog_term_t;

typedef struct isog_args_t {         // mirrors IsoArgs
    md_gpu_float4 origin_h;
    uint32_t num_shells, num_ao, num_fn, fn_stride;
    uint32_t num_fields, level, brick_level, cap;
    uint32_t hash_mask, mode, count, grid_mul;
    uint32_t grid_div;
    float    step_max, lambda, tol;
    uint32_t canon, _pad0;
    uint32_t dform, num_dmat;
    md_gpu_addr_t shells, prims, crit, ao_nrm, ao_ijk, C, terms, fields;
    md_gpu_addr_t nodes_in, nodes_out, count_in, count_out, grid_out;
    md_gpu_addr_t bricks, brick_slot, hash_keys, hash_vals, corners, cell_vert;
    md_gpu_addr_t vpos, vpos_out, vnrm, vval, vinfo, quads, counters, bucket, remap, vdone;
    md_gpu_addr_t out_pos, out_nrm, out_idx;
    md_gpu_addr_t fn_tail;
    md_gpu_addr_t plist;
    md_gpu_addr_t wts;
    md_gpu_addr_t dmat, fdmat, drs, dprop;
    md_gpu_addr_t scan;
    md_gpu_addr_t bstate;
    uint32_t group_base, grid_max;
    uint32_t dsparse, dm_nt, dm_pool, _pad2;
    md_gpu_addr_t dm_tao, dm_map, dm_rows;
    md_gpu_addr_t pts, pval;
    uint32_t pts_stride, pts_field, _pad3, _pad4;
} isog_args_t;

enum {
    ISOG_K_GRID, ISOG_K_SCAN, ISOG_K_OCTREE, ISOG_K_OCTREE_SMALL, ISOG_K_ODFORM, ISOG_K_OEMIT, ISOG_K_BSLOTS, ISOG_K_HASH,
    ISOG_K_POINTS, ISOG_K_CORNERS_FILL, ISOG_K_CLIST_F1V4,
    ISOG_K_GCORN_F16, ISOG_K_GCORN_F32, ISOG_K_GCORN_F64, ISOG_K_GEVAL_F16, ISOG_K_GEVAL_F32, ISOG_K_GEVAL_F64,
    // the wide GEMM for canonical sets of more than 64 functions, by width: 128, 192, 256 columns
    ISOG_K_GCORN_W128, ISOG_K_GCORN_W192, ISOG_K_GCORN_Y256, ISOG_K_GEVAL_W128, ISOG_K_GEVAL_W192, ISOG_K_GEVAL_Y256,
    ISOG_K_DCORN, ISOG_K_DEVAL,
    ISOG_K_FSEED, ISOG_K_FFACE, ISOG_K_FFINISH, ISOG_K_FATOMS,
    ISOG_K_CELLS, ISOG_K_QUADS, ISOG_K_RELAX, ISOG_K_NEWTON, ISOG_K_EVAL1, ISOG_K_PROJ1,
    ISOG_K_BCOUNT, ISOG_K_BSCAN, ISOG_K_BVERTS, ISOG_K_BTRIS,
    ISOG_K_DROWS,
    ISOG_K_PLOAD, ISOG_K_PSTORE,
    ISOG_K_COUNT
};

static md_gpu_kernel_t isog_k[ISOG_K_COUNT];
static const char*     isog_kname[ISOG_K_COUNT];   // their labels

static md_gpu_device_t isog_device;
static bool             isog_failed;
// Held while the kernels are made: a set is made on whatever thread its caller likes (md_gto_iso_gpu_set_create
// makes them, so that the first extraction's thread need not), and two may be made at once
static md_mutex_t       isog_kernel_mutex;
static bool             isog_kernel_mutex_ok;

// The choices, from the measurements of branch gto-iso-gpu-opt (r1 to r31; the reports in its
// history), fixed where one form won everywhere it was measured:
//   - the bricks of the octree's last tested level's children untested: the corners decide, for less than
//     the last level's bound costs (on the RTX 3080 orbitals 0.2-1 ms faster, flooded densities 2.5 ms);
//   - the octree tests levels from 2 on: the root box holds the surfaces, levels 0 and 1 never dropped a
//     node, and their few nodes reaching every AO took a group each in turn (2.9 ms on the M4 Pro);
//   - a brick's high-face corners from the neighbour that owns them; the corners from a point list in the
//     bricks' (Morton) order, 256 points per group;
//   - one function: four corner points per thread (iso_corners_list_f1v4) and the Newton projection in
//     one launch (iso_project_list_f1); more: the tiled GEMM (iso_gcorners_f*, iso_geval_f*), canonical
//     sets of more than 64 functions in one pass of 128, 192 (iso_gcorners_w*: two groups to an RTX 3080
//     SM) or 256 columns (iso_gcorners_y256: the next chunk's rows of C ahead);
//   - the vertices relaxed on the lattice's surface first, then projected, converged within 1e-2 h, a step
//     predicted far below that ending the vertex.
// What stays a switch (md_gto_iso_gpu_switches_set, or GTO_ISO_* in the environment at initialization):
//   GTO_ISO_EVAL_EPS=10|8|6   the AO cutoff 1e-N of the evaluation (and of the octree's culling, whose
//                             bound adds the dropped part): the volume path uses 1e-6 (default 6)
static int isog_eval_eps_exp = 6;
//   GTO_ISO_DFORM=1           sets of densities: corners and vertices as phi^T D phi over the AOs reaching each
//                             group, D from the fields' terms (iso_dcorners, iso_deval). Sets of more functions
//                             than the factored kernels hold (1024) use it regardless; else it is slower than the
//                             factored GEMM (C60 density) (default 0)
static int isog_dform = 0;
//   GTO_ISO_FLOOD=0|1|2       APPROXIMATE: the corners from seeds and a flood instead of every brick the octree
//                             keeps. The bricks' corners are evaluated first; one whose corners straddle an iso
//                             value is a seed, evaluated in full; then, round by round, a brick across a face of
//                             an evaluated one with a sign change on that face. Only bricks with surface cells
//                             are evaluated, but a surface piece no seed reaches (smaller than about 4 cells, or
//                             between the brick corners) is missed. The bricks holding an atom are seeds too, so
//                             a total density's pieces (each holds a nucleus, its maxima are there) are reached
//                             however small. 1: densities whose maxima are at the atoms (SCF densities, total or
//                             of one spin: md_gto_iso_field_set_maxima_at_atoms) with every weight >= 0, exact on
//                             every case measured (iso 0.0025 to 3); 2: every field (to see what it misses:
//                             orbitals' small lobes, attachment and detachment densities' pieces between the
//                             atoms at high iso values); 0 off (default 1)
static int isog_flood = 1;
//   GTO_ISO_FLOOD_ATOMS=0     TEST: the flood without the atoms' bricks as seeds (default 1)
static int isog_flood_atoms = 1;
//   GTO_ISO_CHUNKS=n          TEST: the corner and vertex launches in n chunks, each followed by a progress marker
//                             (default 0: as many as the octree's time suggests, 1 when it was quick)
static int isog_chunks = 0;
//   GTO_ISO_DSPARSE=0         TEST: the D form's matrices dense, made on the host (num_ao^2 floats each: 371 MB for a
//                             protein of 9636 AOs, and its product C W C^T seconds of the CPU), rather than in tiles
//                             of the AO pairs that meet (md_gto_gpu_dmat), built on the device (default 1)
static int isog_dsparse = 1;
// The pattern's margin (md_gto_gpu_dmat_desc_t): the AO lists and the pattern take their radii from different
// bounds; a pair the lists hold and the pattern not is zero in the matrix (both AOs below the cutoff wherever one is)
#define ISOG_DMAT_MARGIN 1.0

#define ISOG_NEWTON_TOL_EXP 2       // converged within 1e-2 h: the facets' own deviation is about 3e-2 h rms
#define ISOG_OCTREE_START   2       // the first level the octree tests

// The kernels are made on first use, so initializing costs nothing until a set is made
void md_gto_iso_gpu_initialize(md_gpu_device_t device) {
    isog_device = device;
    isog_failed = false;
    if (!isog_kernel_mutex_ok) isog_kernel_mutex_ok = md_mutex_init(&isog_kernel_mutex);
    const char* e = getenv("GTO_ISO_EVAL_EPS");
    if (e && e[0]) isog_eval_eps_exp = atoi(e);
    if (isog_eval_eps_exp < 4 || isog_eval_eps_exp > 12) isog_eval_eps_exp = 6;
    e = getenv("GTO_ISO_DFORM");
    if (e && e[0]) isog_dform = atoi(e);
    e = getenv("GTO_ISO_FLOOD");
    if (e && e[0]) isog_flood = atoi(e);
    if (isog_flood < 0 || isog_flood > 2) isog_flood = 1;
    e = getenv("GTO_ISO_FLOOD_ATOMS");
    if (e && e[0]) isog_flood_atoms = atoi(e);
    e = getenv("GTO_ISO_CHUNKS");
    if (e && e[0]) isog_chunks = atoi(e);
    e = getenv("GTO_ISO_DSPARSE");
    if (e && e[0]) isog_dsparse = atoi(e);
}

// Every kernel is made at once. One a device cannot make (more groupshared memory than it allows, say)
// fails the extraction only where it is chosen (the wide GEMM falls back to tiles of 64); the kernels
// every extraction needs fail it always.
static bool isog_ready;

static bool isog_kernel_core(int i) {
    return i == ISOG_K_GRID || i == ISOG_K_SCAN || i == ISOG_K_OEMIT || i == ISOG_K_BSLOTS || i == ISOG_K_HASH ||
           i == ISOG_K_POINTS || i == ISOG_K_CORNERS_FILL || i == ISOG_K_CELLS || i == ISOG_K_QUADS || i == ISOG_K_RELAX ||
           i == ISOG_K_NEWTON || i == ISOG_K_BCOUNT || i == ISOG_K_BSCAN || i == ISOG_K_BVERTS || i == ISOG_K_BTRIS;
}

// Kernel k, or false (with a message) when this device could not make it
static bool isog_have(int k) {
    if (k >= 0 && k < ISOG_K_COUNT && isog_k[k]) return true;
    MD_LOG_ERROR("md_gto_iso_gpu_extract: kernel %s is not available on this device", (k >= 0 && k < ISOG_K_COUNT && isog_kname[k]) ? isog_kname[k] : "?");
    return false;
}

static bool isog_kernels_make(void);

static bool isog_kernels_ready(void) {
    if (isog_kernel_mutex_ok) md_mutex_lock(&isog_kernel_mutex);
    const bool ok = isog_kernels_make();
    if (isog_kernel_mutex_ok) md_mutex_unlock(&isog_kernel_mutex);
    return ok;
}

static bool isog_kernels_make(void) {
    if (isog_ready) return true;
    if (isog_failed) return false;
    if (!isog_device) { MD_LOG_ERROR("md_gto_iso_gpu_extract: md_gto_iso_gpu_initialize not called"); return false; }
    const md_gpu_kernel_desc_t d[ISOG_K_COUNT] = {
        md_shader_gto_iso_gpu_iso_grid_kernel(),
        md_shader_gto_iso_gpu_iso_scan_kernel(),
        md_shader_gto_iso_gpu_iso_octree_kernel(),
        md_shader_gto_iso_gpu_iso_octree_small_kernel(),
        md_shader_gto_iso_gpu_iso_octree_dform_kernel(),
        md_shader_gto_iso_gpu_iso_octree_emit_kernel(),
        md_shader_gto_iso_gpu_iso_brick_slots_kernel(),
        md_shader_gto_iso_gpu_iso_hash_kernel(),
        md_shader_gto_iso_gpu_iso_points_kernel(),
        md_shader_gto_iso_gpu_iso_corners_fill_kernel(),
        md_shader_gto_iso_gpu_iso_corners_list_f1v4_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_f16_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_f32_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_f64_kernel(),
        md_shader_gto_iso_gpu_iso_geval_f16_kernel(),
        md_shader_gto_iso_gpu_iso_geval_f32_kernel(),
        md_shader_gto_iso_gpu_iso_geval_f64_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_w128_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_w192_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_y256_kernel(),
        md_shader_gto_iso_gpu_iso_geval_w128_kernel(),
        md_shader_gto_iso_gpu_iso_geval_w192_kernel(),
        md_shader_gto_iso_gpu_iso_geval_y256_kernel(),
        md_shader_gto_iso_gpu_iso_dcorners_kernel(),
        md_shader_gto_iso_gpu_iso_deval_kernel(),
        md_shader_gto_iso_gpu_iso_flood_seed_kernel(),
        md_shader_gto_iso_gpu_iso_flood_face_kernel(),
        md_shader_gto_iso_gpu_iso_flood_finish_kernel(),
        md_shader_gto_iso_gpu_iso_flood_atoms_kernel(),
        md_shader_gto_iso_gpu_iso_cells_kernel(),
        md_shader_gto_iso_gpu_iso_quads_kernel(),
        md_shader_gto_iso_gpu_iso_relax_kernel(),
        md_shader_gto_iso_gpu_iso_newton_kernel(),
        md_shader_gto_iso_gpu_iso_eval_orbital_kernel(),
        md_shader_gto_iso_gpu_iso_project_list_f1_kernel(),
        md_shader_gto_iso_gpu_iso_bucket_count_kernel(),
        md_shader_gto_iso_gpu_iso_bucket_scan_kernel(),
        md_shader_gto_iso_gpu_iso_bucket_vertices_kernel(),
        md_shader_gto_iso_gpu_iso_bucket_triangles_kernel(),
        md_shader_gto_iso_gpu_iso_dmat_rows_kernel(),
        md_shader_gto_iso_gpu_iso_points_load_kernel(),
        md_shader_gto_iso_gpu_iso_points_store_kernel(),
    };
    for (int i = 0; i < ISOG_K_COUNT; ++i) {
        if (!isog_k[i]) isog_k[i] = md_gpu_kernel_create(isog_device, &d[i]);
        isog_kname[i] = d[i].label;
        if (!isog_k[i]) {
            MD_LOG_ERROR("md_gto_iso_gpu: kernel %s: %s", d[i].label ? d[i].label : "?", md_gpu_last_error() ? md_gpu_last_error() : "?");
            if (isog_kernel_core(i)) isog_failed = true;
        }
    }
    if (isog_failed) return false;
    isog_ready = true;
    return true;
}

void md_gto_iso_gpu_shutdown(void) {
    for (int i = 0; i < ISOG_K_COUNT; ++i) {
        if (isog_k[i]) { md_gpu_kernel_destroy(isog_k[i]); isog_k[i] = NULL; }
    }
    isog_device = NULL;
    isog_failed = false;
    isog_ready = false;
    if (isog_kernel_mutex_ok) md_mutex_destroy(&isog_kernel_mutex);
    isog_kernel_mutex_ok = false;
}

// ---------------------------------------------------------------------------------------------
// The field set
// ---------------------------------------------------------------------------------------------

#define ISOG_FN_CAP 1024     // O_FN_CAP in the shader

struct md_gto_iso_gpu_set_t {
    md_gpu_mem_t shells, prims, crit, ao_nrm, ao_ijk, C, terms, fn_tail, wts, dmat, fdmat, drs, dprop;
    md_gpu_mem_t dm_tao, dm_map, dm_rows;   // dsparse: the pattern the matrices are in (md_gto_gpu_dmat_view)
    uint32_t     dsparse, dm_nt, dm_pool;   // dm_pool: floats per matrix
    uint32_t     num_shells, num_ao, num_fn, fn_stride, num_fields;
    bool         orbitals_only;     // every field one function, one linear term: iso_eval_orbital
    bool         canon;             // every term w a, w a^2 or w a b (a even, b = a + 1): wts holds them per column
    uint32_t     num_dmat;          // > 0: every field a density, as phi^T D phi: dmat (one per distinct field), fdmat
    bool         factorless;        // more functions than the factored kernels hold: the D form only
    double       d_sraw[MD_GTO_ISO_GPU_MAX_FIELDS];    // per matrix: sum of |D| (without the AO normalization)
    float        d_slack[MD_GTO_ISO_GPU_MAX_FIELDS];   // its field's slack
    bool         d_psd[MD_GTO_ISO_GPU_MAX_FIELDS];     // positive semidefinite (every term w a^2, w >= 0)
    isog_field_t fields[MD_GTO_ISO_GPU_MAX_FIELDS];    // iso filled in per extraction
    iso_field_t  host[MD_GTO_ISO_GPU_MAX_FIELDS];         // shallow copies for the root box (their tails): the fields outlive the set
    const md_gto_iso_field_t* src[MD_GTO_ISO_GPU_MAX_FIELDS];
    struct md_allocator_i* alloc;
    isog_shell_t* prof_shells;      // BENCHMARK: the shells on the host, for the profile's work counts
    // the evaluation cutoff the shells' radii (c.w) and fn_tail are for, and what re-derives them
    int           eps_exp;
    double*       sh_nmax;          // per shell: the largest AO normalization
    double*       sh_amin;          // per shell: the smallest exponent
    md_gpu_float2* host_prims;      // alpha, contraction coefficient
    double*       colsum;           // per function column: sum of |C|
    double*       sh_cmax;          // per shell: the largest |C| of its AOs over the functions (coefficient screening)
    int           screen;           // the radii are coefficient-screened (sets of orbitals)
};

static inline float isog_f32_up(double x);

// Every AO of shell s below eps beyond the returned radius (squared)
static double isog_shell_r2(const md_gto_iso_gpu_set_t* set, uint32_t s, double eps) {
    const isog_shell_t* sh = &set->prof_shells[s];
    const uint32_t l = sh->info.x, nprim = sh->info.z, poff = sh->info.w;
    double r = sqrt(l / (2.0 * set->sh_amin[s]));
    for (;; r += 0.05) {
        double t = 0.0;
        for (uint32_t q = 0; q < nprim; ++q) t += fabs((double)set->host_prims[poff + q].y) * exp(-(double)set->host_prims[poff + q].x * r * r);
        if (set->sh_nmax[s] * pow(r, (double)l) * t < eps || r > 1.0e3) break;
    }
    return r * r;
}

// The matrices' tail constants for the cutoff eps (what the AOs a node drops, each below eps over it,
// can add; doubled for the float distances that decide it), uploaded
static void isog_set_dprop(md_gpu_stream_t stream, md_gto_iso_gpu_set_t* set, double eps) {
    if (!set->num_dmat || !set->dprop.gpu) return;
    md_gpu_float4 pr[MD_GTO_ISO_GPU_MAX_FIELDS];
    for (uint32_t u = 0; u < set->num_dmat; ++u) {
        pr[u].x = isog_f32_up(4.0 * eps * eps * set->d_sraw[u] * (1.0 + 1.0e-6));
        pr[u].y = isog_f32_up(4.0 * eps * (1.0 + 1.0e-6));
        pr[u].z = set->d_psd[u] ? 1.0f : 0.0f;
        pr[u].w = set->d_slack[u];
    }
    md_gpu_upload(stream, set->dprop.gpu, pr, sizeof(md_gpu_float4) * set->num_dmat);
}

// The shells' radii and the functions' tails for the cutoff 1e-eps_exp, uploaded
static void isog_set_eps(md_gpu_stream_t stream, md_gto_iso_gpu_set_t* set, int eps_exp, int screen) {
    const double eps = pow(10.0, -(double)eps_exp);
    screen = screen && set->orbitals_only && set->sh_cmax;
    for (uint32_t s = 0; s < set->num_shells; ++s) {
        if (!screen) {
            set->prof_shells[s].c.w = isog_f32_up(isog_shell_r2(set, s, eps));
        } else if (set->sh_cmax[s] > 0.0) {
            // coefficient screening: |C| |AO| below eps beyond the radius
            set->prof_shells[s].c.w = isog_f32_up(isog_shell_r2(set, s, eps / set->sh_cmax[s]));
        } else {
            set->prof_shells[s].c.w = -1.0f;    // no coefficient: never listed
        }
    }
    md_gpu_upload(stream, set->shells.gpu, set->prof_shells, sizeof(isog_shell_t) * set->num_shells);
    // per function: what the shells that do not reach a box can add over it, every AO of theirs being
    // below eps there (screened: every coefficient times its AO), twice that for the float distances
    // that decide it
    float* tail = (float*)calloc(set->fn_stride + 1, sizeof(float));
    for (uint32_t j = 0; j < set->num_fn; ++j) {
        const double t = screen ? 2.0 * eps * (double)set->num_ao : 2.0 * eps * set->colsum[j];
        tail[j] = isog_f32_up(t * (1.0 + 1.0e-6));
    }
    md_gpu_upload(stream, set->fn_tail.gpu, tail, sizeof(float) * set->fn_stride);
    free(tail);
    isog_set_dprop(stream, set, eps);
    set->eps_exp = eps_exp;
    set->screen = screen;
}

static inline float isog_f32_dn(double x) { float f = (float)x; return (double)f > x ? nextafterf(f, -INFINITY) : f; }
static inline float isog_f32_up(double x) { float f = (float)x; return (double)f < x ? nextafterf(f,  INFINITY) : f; }

// A density's matrix from its terms, D = sum_t w_t sym(a_t b_t^T), into Dd (num_ao^2): every term as
// squares (a square as it is, a product as (w/4) ((a + b)^2 - (a - b)^2)) through the blocked product
// md_gto_density_from_orbitals. An SCF density from its orbitals is num_occ squares, which a plain loop
// over the terms would take num_occ passes over the whole matrix for. False when out of memory.
static bool isog_dmat_from_terms(double* Dd, const iso_density_t* d) {
    const size_t na = (size_t)d->num_ao, nf = (size_t)d->num_fn;
    bool squares = true;
    for (int t = 0; t < d->num_terms && squares; ++t) squares = d->terms[t].a == d->terms[t].b;
    if (squares) {
        // per function the sum of its terms' weights: the functions are the columns of d->C
        double* w = (double*)calloc(nf + 1, sizeof(double));
        if (!w) return false;
        for (int t = 0; t < d->num_terms; ++t) w[d->terms[t].a] += d->terms[t].w;
        const bool ok = md_gto_density_from_orbitals(Dd, d->C, 1, nf, w, nf, na);
        free(w);
        return ok;
    }
    size_t nr = 0;
    for (int t = 0; t < d->num_terms; ++t) nr += d->terms[t].a == d->terms[t].b ? 1 : 2;
    double* rows = (double*)malloc(sizeof(double) * (nr * na + 1));
    double* w = (double*)malloc(sizeof(double) * (nr + 1));
    bool ok = rows && w;
    if (ok) {
        size_t r = 0;
        for (int t = 0; t < d->num_terms; ++t) {
            const iso_term_t* tm = &d->terms[t];
            double* p = rows + r * na;
            if (tm->a == tm->b) {
                for (size_t i = 0; i < na; ++i) p[i] = d->C[i * nf + (size_t)tm->a];
                w[r++] = tm->w;
            } else {
                double* q = p + na;
                for (size_t i = 0; i < na; ++i) {
                    const double ca = d->C[i * nf + (size_t)tm->a], cb = d->C[i * nf + (size_t)tm->b];
                    p[i] = ca + cb;
                    q[i] = ca - cb;
                }
                w[r++] = 0.25 * tm->w;
                w[r++] = -0.25 * tm->w;
            }
        }
        ok = md_gto_density_from_orbitals(Dd, rows, na, 1, w, nr, na);
    }
    free(w);
    free(rows);
    return ok;
}

// The D form's matrices in tiles of the AO pairs that meet (md_gto_gpu_dmat), one pattern for all, built on the device
// from each density's functions, or uploaded from the matrix it was made of: set->dmat holds num_dmat pools of
// set->dm_pool floats, set->dm_tao and set->dm_map the pattern. Per field its matrix into fdm, and the matrices'
// properties (d_psd, d_slack; d_sraw with the row sums, isog_dmat_rows). False, with nothing kept, where md_gto
// cannot (md_gto_gpu_initialize not called, out of memory): the matrices are then made dense.
static bool isog_dmat_tiles(md_gto_iso_gpu_set_t* set, md_gpu_stream_t stream, const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
                            const md_gto_iso_field_t* const* fields, uint32_t num_fields, const int* first_of, uint32_t* fdm, uint32_t num_dmat) {
    // for the cutoff the extraction's AO lists are for (isog_set_eps), plus a margin for their different radii
    const md_gto_gpu_dmat_desc_t desc = {
        .basis = basis, .atom_xyz = atom_xyz, .atom_xyz_stride = atom_xyz_stride,
        .cutoff = pow(10.0, -(double)isog_eval_eps_exp), .margin = ISOG_DMAT_MARGIN,
    };
    const md_tick_t t0 = md_tick_now();
    md_gto_gpu_dmat_t dm = md_gto_gpu_dmat_create(stream, &desc);
    if (!dm) return false;
    const double t_pattern = md_tick_to_milliseconds(md_tick_now() - t0);
    md_gto_gpu_dmat_view_t v;
    md_gto_gpu_dmat_view(&v, dm);
    const size_t nao = v.num_ao;
    const size_t pool = 256 * (size_t)MAX(1, v.num_tiles);
    bool ok = pool * num_dmat < (size_t)UINT32_MAX;
    if (ok) {
        set->dmat   = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(float) * pool * num_dmat);
        set->dm_tao = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(uint32_t) * MAX(1, nao));
        set->dm_map = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(uint32_t) * (size_t)v.nt * v.nt);
        ok = set->dmat.gpu && set->dm_tao.gpu && set->dm_map.gpu;
    }
    ok = ok && md_gpu_copy(stream, set->dm_tao.gpu, v.tao, sizeof(uint32_t) * nao)
            && md_gpu_copy(stream, set->dm_map.gpu, v.map, sizeof(uint32_t) * (size_t)v.nt * v.nt);
    uint32_t u = 0;
    for (uint32_t k = 0; k < num_fields && ok; ++k) {
        if (first_of[k] != (int)k) { fdm[k] = fdm[first_of[k]]; continue; }
        const md_gto_iso_field_t* fk = fields[k];
        const iso_density_t* d = fk->f.dens;
        bool squares = true;
        for (int t = 0; t < d->num_terms && squares; ++t) squares = d->terms[t].a == d->terms[t].b;
        if (fk->dmat && fk->dmat_dim == nao) {
            ok = md_gto_gpu_dmat_set_matrix(stream, dm, fk->dmat);              // the matrix the field was made from
        } else if (squares) {
            // per function the sum of its terms' weights: the functions are the columns of d->C
            double* w = (double*)calloc((size_t)d->num_fn + 1, sizeof(double));
            ok = w != NULL;
            for (int t = 0; ok && t < d->num_terms; ++t) w[d->terms[t].a] += d->terms[t].w;
            ok = ok && md_gto_gpu_dmat_build_columns(stream, dm, d->C, (size_t)d->num_fn, w, (size_t)d->num_fn);
            free(w);
        } else {
            // a product as (w/4) ((a + b)^2 - (a - b)^2), as isog_dmat_from_terms does
            size_t nr = 0;
            for (int t = 0; t < d->num_terms; ++t) nr += d->terms[t].a == d->terms[t].b ? 1 : 2;
            const size_t nf = (size_t)d->num_fn;
            double* rows = (double*)malloc(sizeof(double) * (nr * nao + 1));
            double* w = (double*)malloc(sizeof(double) * (nr + 1));
            ok = rows && w;
            size_t r = 0;
            for (int t = 0; ok && t < d->num_terms; ++t) {
                const iso_term_t* tm = &d->terms[t];
                double* p = rows + r * nao;
                if (tm->a == tm->b) {
                    for (size_t i = 0; i < nao; ++i) p[i] = d->C[i * nf + (size_t)tm->a];
                    w[r++] = tm->w;
                } else {
                    double* q = p + nao;
                    for (size_t i = 0; i < nao; ++i) {
                        const double ca = d->C[i * nf + (size_t)tm->a], cb = d->C[i * nf + (size_t)tm->b];
                        p[i] = ca + cb;
                        q[i] = ca - cb;
                    }
                    w[r++] = 0.25 * tm->w;
                    w[r++] = -0.25 * tm->w;
                }
            }
            ok = ok && md_gto_gpu_dmat_build(stream, dm, rows, w, nr);
            free(rows);
            free(w);
        }
        ok = ok && md_gpu_copy(stream, set->dmat.gpu + sizeof(float) * pool * u, v.tiles, sizeof(float) * pool);
        bool psd = true;
        for (int t = 0; t < d->num_terms && psd; ++t) psd = d->terms[t].a == d->terms[t].b && d->terms[t].w >= 0.0;
        set->d_psd[u] = psd;
        set->d_slack[u] = isog_f32_up(d->slack);
        fdm[k] = u++;
    }
    md_gto_gpu_dmat_destroy(stream, dm);    // stream-ordered: after the copies
    if (!ok || u != num_dmat) {
        md_gpu_free(stream, set->dmat.gpu); md_gpu_free(stream, set->dm_tao.gpu); md_gpu_free(stream, set->dm_map.gpu);
        set->dmat = set->dm_tao = set->dm_map = (md_gpu_mem_t){0};
        MD_LOG_DEBUG("md_gto_iso_gpu_set_create: the matrices in tiles could not be made; dense");
        return false;
    }
    set->dsparse = 1;
    set->dm_nt = v.nt;
    set->dm_pool = (uint32_t)pool;
    MD_LOG_DEBUG("md_gto_iso_gpu_set_create: %u matrices in tiles, %llu of %llu tile pairs (%.1f MB each); the pattern %.1f ms, the builds' host part %.1f ms",
                 num_dmat, (unsigned long long)v.num_tiles, (unsigned long long)v.nt * (v.nt + 1) / 2, sizeof(float) * pool / 1.0e6, t_pattern,
                 md_tick_to_milliseconds(md_tick_now() - t0) - t_pattern);
    return true;
}

// The tiles' row sums (iso_dmat_rows): the octree's tails (drs), and per matrix the sum of |D| (d_sraw), read back.
// After the AO normalizations are on the device.
static bool isog_dmat_rows(md_gto_iso_gpu_set_t* set, md_gpu_stream_t stream) {
    if (!isog_have(ISOG_K_DROWS)) return false;
    const uint32_t n = set->num_ao, nm = set->num_dmat;
    set->drs = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(float) * (size_t)nm * n);
    set->dm_rows = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(float) * (size_t)nm * n);
    md_gpu_mem_t rb = md_gpu_malloc(stream, MD_GPU_MEM_HOST_READ, sizeof(float) * (size_t)nm * n);
    if (!set->drs.gpu || !set->dm_rows.gpu || !rb.cpu) {
        MD_LOG_ERROR("md_gto_iso_gpu_set_create: out of device memory for the matrices' row sums");
        md_gpu_free(stream, rb.gpu);
        return false;
    }
    isog_args_t a;
    memset(&a, 0, sizeof(a));
    a.num_ao = n;
    a.num_dmat = nm;
    a.dsparse = 1;
    a.dm_nt = set->dm_nt;
    a.dm_pool = set->dm_pool;
    a.dmat = set->dmat.gpu;
    a.dm_tao = set->dm_tao.gpu;
    a.dm_map = set->dm_map.gpu;
    a.dm_rows = set->dm_rows.gpu;
    a.ao_nrm = set->ao_nrm.gpu;
    a.drs = set->drs.gpu;
    bool ok = md_gpu_launch(stream, isog_k[ISOG_K_DROWS], md_gpu_grid((n + 63u) / 64u, nm, 1), &a, sizeof(a));
    ok = ok && md_gpu_copy(stream, rb.gpu, set->dm_rows.gpu, sizeof(float) * (size_t)nm * n);
    md_gpu_stream_sync(stream);
    for (uint32_t u = 0; ok && u < nm; ++u) {
        const float* r = (const float*)rb.cpu + (size_t)u * n;
        double t = 0.0;
        for (uint32_t i = 0; i < n; ++i) t += (double)r[i];
        set->d_sraw[u] = t * 1.001;     // the float sums' rounding, as drs
    }
    md_gpu_free(stream, rb.gpu);
    return ok;
}

md_gto_iso_gpu_set_t* md_gto_iso_gpu_set_create(md_gpu_stream_t stream, const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
                                                const md_gto_iso_field_t* const* fields, uint32_t num_fields, struct md_allocator_i* alloc) {
    if (!stream || !basis || !atom_xyz || !fields || !alloc || num_fields == 0 || num_fields > MD_GTO_ISO_GPU_MAX_FIELDS) {
        MD_LOG_ERROR("md_gto_iso_gpu_set_create: invalid arguments (1 to %d fields)", MD_GTO_ISO_GPU_MAX_FIELDS);
        return NULL;
    }
    // The kernels with the set, on its thread: the first extraction then costs its own thread nothing for them
    // (seconds on a CPU device, a fraction of one on most GPUs)
    if (!isog_kernels_ready()) return NULL;
    size_t num_atoms;
    float* xyz = iso_pack_atoms(basis, atom_xyz, atom_xyz_stride, &num_atoms);
    iso_density_t B = {0};      // the basis as contracted shells, and their radial functions
    bool ok = iso_density_ao_shells(&B, basis, xyz, NULL);
    free(xyz);
    if (!ok) { free(B.ao_shells); MD_LOG_ERROR("md_gto_iso_gpu_set_create: unsupported basis"); return NULL; }
    iso_radial_init(&B.radial, B.ao_shells, B.num_ao_shells);
    const uint32_t nao = (uint32_t)B.num_ao;

    // functions: every field's, each field starting at an even index (pairs never straddle a tile).
    // A field given more than once (one density at several iso values) shares its functions.
    int first_of[MD_GTO_ISO_GPU_MAX_FIELDS];
    for (uint32_t k = 0; k < num_fields; ++k) {
        first_of[k] = (int)k;
        for (uint32_t j = 0; j < k; ++j) if (fields[j] == fields[k]) { first_of[k] = (int)j; break; }
    }
    // a set of orbitals has no pair terms, so its functions need no alignment
    bool orbitals_only = true;
    for (uint32_t k = 0; k < num_fields; ++k) orbitals_only &= fields[k] && fields[k]->f.dens == NULL;
    const uint32_t align = orbitals_only ? 0u : 1u;
    uint32_t nf = 0;
    for (uint32_t k = 0; k < num_fields; ++k) {
        const md_gto_iso_field_t* fk = fields[k];
        if (!fk) { ok = false; break; }
        if (first_of[k] != (int)k) continue;
        const uint32_t n = fk->f.dens ? (uint32_t)fk->f.dens->num_fn : 1u;
        if ((fk->f.dens ? (uint32_t)fk->f.dens->num_ao : (uint32_t)fk->num_coeff) != nao) { ok = false; break; }
        nf = (nf + align) & ~align;
        nf += n;
    }
    // the D form (phi^T D phi): every field a density of products (no linear terms), the indices of the
    // local lists in 16 bits. With more functions than the factored kernels hold, the set has no
    // functions at all and only the D form runs.
    bool dform_ok = ok && B.num_ao_shells <= 65535 && nao <= 65535;
    for (uint32_t k = 0; k < num_fields && dform_ok; ++k) {
        const iso_density_t* d = fields[k]->f.dens;
        if (!d) { dform_ok = false; break; }
        for (int t = 0; t < d->num_terms; ++t) if (d->terms[t].b < 0) dform_ok = false;
    }
    // TEST: GTO_ISO_FN_CAP lowers the limit (to exercise the factorless sets on small systems)
    uint32_t fn_cap = ISOG_FN_CAP;
    const char* ce = getenv("GTO_ISO_FN_CAP");
    if (ce && atoi(ce) > 0 && (uint32_t)atoi(ce) < fn_cap) fn_cap = (uint32_t)atoi(ce);
    const bool factorless = nf > fn_cap && dform_ok;
    if (!ok || (nf > fn_cap && !factorless)) {
        md_array_free(B.radial.shells, md_get_heap_allocator()); md_array_free(B.radial.crit, md_get_heap_allocator()); free(B.ao_shells);
        MD_LOG_ERROR("md_gto_iso_gpu_set_create: fields of a different basis, or more than %d functions in all (and not all densities)", ISOG_FN_CAP);
        return NULL;
    }
    if (factorless) nf = 0;
    const uint32_t stride = nf > 0 ? ((nf + 63u) / 64u) * 64u : 64u;
    float* C = (float*)calloc((size_t)nao * stride + 1, sizeof(float));
    isog_term_t* terms = NULL;
    md_allocator_i* heap = md_get_heap_allocator();

    md_gto_iso_gpu_set_t* set = (md_gto_iso_gpu_set_t*)md_alloc(alloc, sizeof(md_gto_iso_gpu_set_t));
    memset(set, 0, sizeof(*set));
    set->alloc = alloc;
    uint32_t f0 = 0;
    for (uint32_t k = 0; k < num_fields && ok; ++k) {
        const md_gto_iso_field_t* fk = fields[k];
        set->src[k] = fk;
        set->host[k] = fk->f;
        if (first_of[k] != (int)k) { set->fields[k] = set->fields[first_of[k]]; continue; }
        f0 = (f0 + align) & ~align;
        set->fields[k].term_off = (uint32_t)md_array_size(terms);
        if (fk->f.dens && factorless) {
            set->fields[k].slack = isog_f32_up(fk->f.dens->slack);    // no functions, no terms
        } else if (fk->f.dens) {
            const iso_density_t* d = fk->f.dens;
            for (uint32_t i = 0; i < nao; ++i) for (int j = 0; j < d->num_fn; ++j) C[(size_t)i * stride + f0 + j] = (float)d->C[(size_t)i * d->num_fn + j];
            for (int t = 0; t < d->num_terms; ++t) {
                const iso_term_t* tm = &d->terms[t];
                isog_term_t g = {(float)tm->w, f0 + (uint32_t)tm->a, tm->b < 0 ? 0xFFFFFFFFu : f0 + (uint32_t)tm->b, 0};
                if (tm->b >= 0 && (g.a / 64u) != (g.b / 64u)) ok = false;   // the epilogue sees one tile of 64
                md_array_push(terms, g, heap);
            }
            set->fields[k].slack = isog_f32_up(d->slack);
            f0 += (uint32_t)d->num_fn;
        } else {
            for (uint32_t i = 0; i < nao; ++i) C[(size_t)i * stride + f0] = (float)fk->coeff[i];
            isog_term_t g = {1.0f, f0, 0xFFFFFFFFu, 0};
            md_array_push(terms, g, heap);
            f0 += 1u;
        }
        set->fields[k].term_count = (uint32_t)md_array_size(terms) - set->fields[k].term_off;
    }

    // the terms per field and column, when each is w a, w a^2 or w a b with a even and b = a + 1 (what
    // orbitals, density matrices and products of pairs give): the register-form epilogues' weights
    bool canon = ok;
    float* W = (float*)calloc((size_t)num_fields * 3 * stride + 1, sizeof(float));
    for (uint32_t k = 0; k < num_fields && canon; ++k) {
        float* wl = W + (size_t)k * 3 * stride;
        float* wd = wl + stride;
        float* wp = wd + stride;
        for (uint32_t t = 0; t < set->fields[k].term_count; ++t) {
            const isog_term_t* g = &terms[set->fields[k].term_off + t];
            uint32_t ta = g->a, tb = g->b;
            if (tb == 0xFFFFFFFFu) { wl[ta] += g->w; continue; }
            if (ta == tb) { wd[ta] += g->w; continue; }
            if (tb < ta) { const uint32_t x = ta; ta = tb; tb = x; }
            if ((ta & 1u) == 0 && tb == ta + 1) wp[ta] += g->w;
            else canon = false;
        }
    }

    // shells, primitives, radial critical intervals (outward rounded), AO tables
    isog_shell_t* sh = (isog_shell_t*)calloc(B.num_ao_shells + 1, sizeof(isog_shell_t));
    double* sh_nmax = (double*)calloc(B.num_ao_shells + 1, sizeof(double));
    double* sh_amin = (double*)calloc(B.num_ao_shells + 1, sizeof(double));
    md_gpu_float2* pr = NULL;
    md_gpu_float4* cr = NULL;
    float*    nrm = (float*)calloc(nao + 1, sizeof(float));
    uint32_t* ijk = (uint32_t*)calloc(nao + 1, sizeof(uint32_t));
    for (int s = 0; s < B.num_ao_shells; ++s) {
        const iso_ao_shell_t* a = &B.ao_shells[s];
        const iso_shell_t* g = &B.radial.shells[s];
        if (a->nprim > 64) ok = false;      // the corner and vertex kernels stage at most 64 per shell
        // point values: every AO of the shell below 1e-10 beyond this radius (the corner values only
        // decide on which side of an iso value a corner lies; the bounds do not use it)
        double nmax = 0.0, amin = DBL_MAX;
        for (int n = 0; n < a->ncart; ++n) nmax = fmax(nmax, fabs(a->nrm[n]));
        for (int q = 0; q < a->nprim; ++q) amin = fmin(amin, a->alpha[q]);
        sh_nmax[s] = nmax;
        sh_amin[s] = amin;
        double r = sqrt(a->l / (2.0 * amin));
        for (;; r += 0.05) {
            double t = 0.0;
            for (int q = 0; q < a->nprim; ++q) t += fabs(a->coef[q]) * exp(-a->alpha[q] * r * r);
            if (nmax * pow(r, a->l) * t < 1.0e-10) break;
        }
        sh[s].c     = (md_gpu_float4){(float)a->c[0], (float)a->c[1], (float)a->c[2], isog_f32_up(r * r)};
        sh[s].info  = (md_gpu_uint4){(uint32_t)a->l, (uint32_t)a->ncart, (uint32_t)a->nprim, (uint32_t)md_array_size(pr)};
        sh[s].info2 = (md_gpu_uint4){(uint32_t)a->ao, (uint32_t)md_array_size(cr), (uint32_t)g->crit_num[0][0], 0};
        for (int q = 0; q < a->nprim; ++q) { md_gpu_float2 p = {(float)a->alpha[q], (float)a->coef[q]}; md_array_push(pr, p, heap); }
        for (int k = 0; k < g->crit_num[0][0]; ++k) {
            const iso_crit_t* c = &B.radial.crit[g->crit_off[0][0] + k];
            md_gpu_float4 cf = {isog_f32_dn(c->lo), isog_f32_up(c->hi), isog_f32_dn(c->g.lo - 1e-7 * fabs(c->g.lo)), isog_f32_up(c->g.hi + 1e-7 * fabs(c->g.hi))};
            md_array_push(cr, cf, heap);
        }
        for (int n = 0; n < a->ncart; ++n) {
            nrm[a->ao + n] = (float)a->nrm[n];
            ijk[a->ao + n] = (uint32_t)a->ijk[n][0] | ((uint32_t)a->ijk[n][1] << 4) | ((uint32_t)a->ijk[n][2] << 8);
        }
    }
    // the densities as matrices (dform): D = sum of the terms' symmetrized products, times the AO
    // normalizations, one per distinct field; when every field is a density of products (no linear
    // terms) and the local lists of the kernels hold every shell
    uint32_t num_dmat = 0;
    float* Dm = NULL;
    float* Rs = NULL;
    uint32_t fdm[MD_GTO_ISO_GPU_MAX_FIELDS] = {0};
    bool dform = dform_ok;
    for (uint32_t k = 0; k < num_fields && dform; ++k) if (first_of[k] == (int)k) num_dmat++;
    // in tiles of the AO pairs that meet, built on the device (isog_dmat_tiles); dense on the host where that fails
    const bool dsparse = dform && isog_dsparse && isog_dmat_tiles(set, stream, basis, atom_xyz, atom_xyz_stride, fields, num_fields, first_of, fdm, num_dmat);
    if (dform && !dsparse) {
        Dm = (float*)calloc((size_t)num_dmat * nao * nao + 1, sizeof(float));
        Rs = (float*)calloc((size_t)num_dmat * nao + 1, sizeof(float));
        double* Dd = (double*)malloc(sizeof(double) * ((size_t)nao * nao + 1));
        uint32_t u = 0;
        for (uint32_t k = 0; k < num_fields && Dm && Rs && Dd; ++k) {
            if (first_of[k] != (int)k) { fdm[k] = fdm[first_of[k]]; continue; }
            const md_gto_iso_field_t* fk = fields[k];
            const iso_density_t* d = fk->f.dens;
            if (fk->dmat && fk->dmat_dim == nao) {
                // the matrix the field was made from
                for (size_t e = 0; e < (size_t)nao * nao; ++e) Dd[e] = fk->dmat[e];
            } else if (!isog_dmat_from_terms(Dd, d)) {
                break;      // out of memory: u stays short of num_dmat
            }
            // with the AO normalizations (the kernels' AO values leave them out); the octree's tails:
            // per AO the sum of |D| over the other AOs (each below eps over a box that drops it), and
            // the sum of |D| without normalizations
            float* out = Dm + (size_t)u * nao * nao;
            double sraw = 0.0;
            for (uint32_t i = 0; i < nao; ++i) {
                double rs = 0.0;
                for (uint32_t j = 0; j < nao; ++j) {
                    const double v = Dd[(size_t)i * nao + j];
                    const float  g = (float)(v * nrm[i] * nrm[j]);
                    out[(size_t)i * nao + j] = g;
                    sraw += fabs(v);
                    rs += fabs((double)g) / fmax(fabs((double)nrm[j]), 1e-300);
                }
                Rs[(size_t)u * nao + i] = isog_f32_up(rs * (1.0 + 1.0e-6));
            }
            bool psd = true;
            for (int t = 0; t < d->num_terms && psd; ++t) psd = d->terms[t].a == d->terms[t].b && d->terms[t].w >= 0.0;
            set->d_sraw[u] = sraw;
            set->d_psd[u] = psd;
            set->d_slack[u] = isog_f32_up(d->slack);
            fdm[k] = u++;
        }
        free(Dd);
        if (!Dm || !Rs || u < num_dmat) dform = false;
    }
    if (!dform) num_dmat = 0;
    set->factorless = factorless;

    set->num_shells = (uint32_t)B.num_ao_shells;
    set->num_ao = nao;
    set->num_fn = nf;
    set->fn_stride = stride;
    set->num_fields = num_fields;
    set->orbitals_only = orbitals_only;
    set->canon = canon;

    // per function column: the sum of |C| (the tails for the shells a box drops: isog_set_eps)
    double* colsum = (double*)calloc(stride + 1, sizeof(double));
    for (uint32_t j = 0; j < nf; ++j) {
        double t = 0.0;
        for (uint32_t i = 0; i < nao; ++i) t += fabs((double)C[(size_t)i * stride + j]);
        colsum[j] = t;
    }
    float* tail = (float*)calloc(stride + 1, sizeof(float));
    for (uint32_t j = 0; j < nf; ++j) tail[j] = isog_f32_up(2.0e-10 * colsum[j] * (1.0 + 1.0e-6));

#define ISOG_UP(mem, ptr, bytes) do { const size_t nb_ = (bytes); set->mem = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, nb_ > 0 ? nb_ : 16); if (set->mem.gpu && nb_ > 0) md_gpu_upload(stream, set->mem.gpu, (ptr), nb_); ok &= set->mem.gpu != 0; } while (0)
    if (ok) {
        ISOG_UP(shells, sh, sizeof(isog_shell_t) * B.num_ao_shells);
        ISOG_UP(prims,  pr, sizeof(md_gpu_float2) * md_array_size(pr));
        ISOG_UP(crit,   cr, sizeof(md_gpu_float4) * md_array_size(cr));
        ISOG_UP(ao_nrm, nrm, sizeof(float) * nao);
        ISOG_UP(ao_ijk, ijk, sizeof(uint32_t) * nao);
        ISOG_UP(C,      C,   sizeof(float) * (size_t)nao * stride);
        ISOG_UP(terms,  terms, sizeof(isog_term_t) * md_array_size(terms));
        ISOG_UP(fn_tail, tail, sizeof(float) * stride);
        ISOG_UP(wts, W, sizeof(float) * (size_t)num_fields * 3 * stride);
        if (num_dmat && !dsparse) {
            ISOG_UP(dmat, Dm, sizeof(float) * (size_t)num_dmat * nao * nao);
            ISOG_UP(drs, Rs, sizeof(float) * (size_t)num_dmat * nao);
        }
        if (num_dmat) {
            ISOG_UP(fdmat, fdm, sizeof(fdm));
            set->dprop = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(md_gpu_float4) * MD_GTO_ISO_GPU_MAX_FIELDS);
            ok &= set->dprop.gpu != 0;
        }
        // the tiles' row sums (drs, and d_sraw), with the AO normalizations above
        set->num_dmat = num_dmat;
        if (num_dmat && dsparse) ok &= isog_dmat_rows(set, stream);
    }
    set->num_dmat = num_dmat;
    if (ok && num_dmat) isog_set_dprop(stream, set, 1.0e-10);    // the radii above are for 1e-10 (isog_set_eps re-derives)
    if (factorless && !num_dmat) ok = false;
    free(Dm);
    free(Rs);
    free(tail);
    free(W);
#undef ISOG_UP
    set->prof_shells = sh;          // BENCHMARK: kept for the profile's work counts
    set->eps_exp = 10;
    set->sh_nmax = sh_nmax;
    set->sh_amin = sh_amin;
    set->colsum = colsum;
    set->host_prims = (md_gpu_float2*)malloc(sizeof(md_gpu_float2) * (md_array_size(pr) + 1));
    memcpy(set->host_prims, pr, sizeof(md_gpu_float2) * md_array_size(pr));
    // per shell, the largest coefficient of its AOs over the set's functions (orbitals: coefficient screening)
    if (orbitals_only) {
        set->sh_cmax = (double*)calloc(B.num_ao_shells + 1, sizeof(double));
        for (int k = 0; k < B.num_ao_shells && set->sh_cmax; ++k) {
            double m = 0.0;
            for (uint32_t n = 0; n < sh[k].info.y; ++n) {
                const size_t ao = sh[k].info2.x + n;
                for (uint32_t j = 0; j < nf; ++j) m = fmax(m, fabs((double)C[ao * stride + j]));
            }
            set->sh_cmax[k] = m;
        }
    }
    free(nrm); free(ijk); free(C);
    md_array_free(pr, heap); md_array_free(cr, heap); md_array_free(terms, heap);
    md_array_free(B.radial.shells, heap); md_array_free(B.radial.crit, heap); free(B.ao_shells);
    if (!ok) {
        MD_LOG_ERROR("md_gto_iso_gpu_set_create: a pair term spans two tiles of 64 functions, a shell of more than 64 primitives, or out of device memory");
        md_gto_iso_gpu_set_destroy(stream, set);
        return NULL;
    }
    return set;
}

void md_gto_iso_gpu_set_destroy(md_gpu_stream_t stream, md_gto_iso_gpu_set_t* set) {
    if (!set) return;
    md_gpu_free(stream, set->shells.gpu); md_gpu_free(stream, set->prims.gpu); md_gpu_free(stream, set->crit.gpu);
    md_gpu_free(stream, set->ao_nrm.gpu); md_gpu_free(stream, set->ao_ijk.gpu); md_gpu_free(stream, set->C.gpu);
    md_gpu_free(stream, set->terms.gpu); md_gpu_free(stream, set->fn_tail.gpu); md_gpu_free(stream, set->wts.gpu);
    md_gpu_free(stream, set->dmat.gpu); md_gpu_free(stream, set->fdmat.gpu);
    md_gpu_free(stream, set->drs.gpu); md_gpu_free(stream, set->dprop.gpu);
    md_gpu_free(stream, set->dm_tao.gpu); md_gpu_free(stream, set->dm_map.gpu); md_gpu_free(stream, set->dm_rows.gpu);
    free(set->prof_shells);
    free(set->sh_nmax); free(set->sh_amin); free(set->colsum); free(set->host_prims); free(set->sh_cmax);
    md_free(set->alloc, set, sizeof(*set));
}

// ---------------------------------------------------------------------------------------------
// Extraction
// ---------------------------------------------------------------------------------------------

// Tests: when set, the extraction copies its bricks, slots and corner values here.
typedef struct isog_debug_t {
    uint32_t  num_bricks, num_slots;
    uint32_t* bricks;       // 4 per brick: x, y, z, sheet mask
    uint32_t* slots;        // per brick
    float*    corners;      // 125 per slot
    double    origin[3];
    int       brick_level;
} isog_debug_t;
static isog_debug_t* isog_debug = NULL;

static bool isog_read(md_gpu_stream_t s, void* dst, md_gpu_addr_t src, size_t bytes) {
    if (!bytes) return true;
    md_gpu_mem_t m = md_gpu_malloc(s, MD_GPU_MEM_HOST_READ, bytes);
    if (!m.gpu) return false;
    md_gpu_copy(s, m.gpu, src, bytes);
    md_gpu_stream_sync(s);
    memcpy(dst, m.cpu, bytes);
    md_gpu_free(s, m.gpu);
    return true;
}

// BENCHMARK: the AOs of the shells that reach the box [lo, hi] (what a kernel lists for it)
static uint32_t isog_prof_box_aos(const md_gto_iso_gpu_set_t* set, const float lo[3], const float hi[3]) {
    uint32_t n = 0;
    for (uint32_t k = 0; k < set->num_shells; ++k) {
        const isog_shell_t* sh = &set->prof_shells[k];
        const float c[3] = {sh->c.x, sh->c.y, sh->c.z};
        float d2 = 0.0f;
        for (int ax = 0; ax < 3; ++ax) {
            const float d = fmaxf(fmaxf(lo[ax] - c[ax], c[ax] - hi[ax]), 0.0f);
            d2 += d * d;
        }
        if (d2 <= sh->c.w) n += sh->info.y;
    }
    return n;
}

// BENCHMARK: the AOs that reach the boxes of consecutive groups of gs points (xyz)
static void isog_prof_groups(const md_gto_iso_gpu_set_t* set, const float* xyz, size_t stride, uint32_t n, uint32_t gs,
                             uint64_t* ao_sum, uint32_t* count) {
    *ao_sum = 0;
    *count = 0;
    for (uint32_t g = 0; g < n; g += gs) {
        float lo[3] = {FLT_MAX, FLT_MAX, FLT_MAX}, hi[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
        for (uint32_t i = g; i < g + gs && i < n; ++i) {
            const float* p = xyz + stride * i;
            for (int ax = 0; ax < 3; ++ax) { lo[ax] = fminf(lo[ax], p[ax]); hi[ax] = fmaxf(hi[ax], p[ax]); }
        }
        *ao_sum += isog_prof_box_aos(set, lo, hi);
        (*count)++;
    }
}

static md_gpu_grid_t isog_groups(uint32_t n) {     // one group per item, split over y
    if (n == 0) return md_gpu_grid(0, 0, 0);
    return md_gpu_grid(n < 65535u ? n : 65535u, (n + 65534u) / 65535u, 1);
}

// The set's cutoff and screening as switched (remade when the switch changed since the set was made)
static void isog_set_current(md_gpu_stream_t s, const md_gto_iso_gpu_set_t* set) {
    const int screen = set->orbitals_only && set->sh_cmax;    // orbitals: a shell dropped where its largest coefficient times its AOs is below eps
    if (set->eps_exp != isog_eval_eps_exp || set->screen != screen) isog_set_eps(s, (md_gto_iso_gpu_set_t*)set, isog_eval_eps_exp, screen);
}

// What every evaluation of the set reads: the basis, the functions (C), the fields' terms and the D form
static void isog_args_set(isog_args_t* a, const md_gto_iso_gpu_set_t* set) {
    a->num_shells = set->num_shells;
    a->num_ao     = set->num_ao;
    a->num_fn     = set->num_fn;
    a->fn_stride  = set->fn_stride;
    a->num_fields = set->num_fields;
    a->shells = set->shells.gpu; a->prims = set->prims.gpu; a->crit = set->crit.gpu;
    a->ao_nrm = set->ao_nrm.gpu; a->ao_ijk = set->ao_ijk.gpu; a->C = set->C.gpu; a->terms = set->terms.gpu;
    a->fn_tail = set->fn_tail.gpu;
    a->wts = set->wts.gpu;
    a->canon = set->canon ? 1u : 0u;
    a->dmat = set->dmat.gpu;
    a->fdmat = set->fdmat.gpu;
    a->num_dmat = set->num_dmat;
    a->dform = ((isog_dform || set->factorless) && set->num_dmat > 0) ? 1u : 0u;
    a->drs = set->drs.gpu;
    a->dsparse = set->dsparse;
    a->dm_nt = set->dm_nt;
    a->dm_pool = set->dm_pool;
    a->dm_tao = set->dm_tao.gpu;
    a->dm_map = set->dm_map.gpu;
    a->dprop = set->dprop.gpu;
}

// The kernel that evaluates value and gradient at a list of vertices (vpos; vinfo.z the field; vdone 0 for the
// ones to evaluate) into vval, and how many vertices a group of it takes (0: a thread each, md_gpu_grid_for).
// 'project': a set of one orbital is evaluated and projected in one launch (iso_project_list_f1), as the
// extraction does; otherwise only evaluated.
//   orbitals:                          a thread per vertex, its field's one column of C (iso_eval_orbital)
//   densities as phi^T D phi (dform):  32 per group (iso_deval)
//   canonical sets of > 64 functions:  the wide GEMM, 16 per group (iso_geval_w128, _w192, _y256)
//   otherwise:                         the tiled GEMM, 32 per group (iso_geval_f16, _f32, _f64)
typedef struct isog_vkernel_t {
    int      k;
    uint32_t per_group;
} isog_vkernel_t;

static isog_vkernel_t isog_vertex_kernel(const md_gto_iso_gpu_set_t* set, bool dform, bool project) {
    if (project && set->orbitals_only && set->num_fn == 1 && set->num_fields == 1) return (isog_vkernel_t){ISOG_K_PROJ1, 128u};
    if (set->orbitals_only) return (isog_vkernel_t){ISOG_K_EVAL1, 0u};
    if (dform) return (isog_vkernel_t){ISOG_K_DEVAL, 32u};
    const uint32_t nf = set->num_fn;
    const int kw = nf <= 128 ? ISOG_K_GEVAL_W128 : nf <= 192 ? ISOG_K_GEVAL_W192 : ISOG_K_GEVAL_Y256;
    if (set->canon && nf > 64u && isog_k[kw]) return (isog_vkernel_t){kw, 16u};
    return (isog_vkernel_t){nf <= 16 ? ISOG_K_GEVAL_F16 : nf <= 32 ? ISOG_K_GEVAL_F32 : ISOG_K_GEVAL_F64, 32u};
}

static void isog_stage(md_gpu_stream_t s, const md_gto_iso_gpu_desc_t* desc, md_gto_iso_gpu_mesh_t* out, int stage, md_tick_t* t) {
    if (!desc->time_stages) return;
    md_gpu_stream_sync(s);
    const md_tick_t now = md_tick_now();
    out->seconds[stage] += md_tick_to_seconds(now - *t);
    *t = now;
}

// ---------------------------------------------------------------------------------------------
// BENCHMARK (temporary, branch gto-iso-gpu-opt): the step by step profile (md_gto_iso_gpu_profile_set)
// ---------------------------------------------------------------------------------------------

static md_gto_iso_gpu_profile_t* isog_prof = NULL;
static md_tick_t isog_prof_tick;

void md_gto_iso_gpu_profile_set(md_gto_iso_gpu_profile_t* profile) { isog_prof = profile; }

void md_gto_iso_gpu_switches(md_gto_iso_gpu_switches_t* out) {
    if (!out) return;
    out->eval_eps_exp = isog_eval_eps_exp;
    out->dform = isog_dform != 0;
    out->flood = isog_flood;
}

void md_gto_iso_gpu_switches_set(const md_gto_iso_gpu_switches_t* in) {
    if (!in) return;
    if (in->eval_eps_exp >= 4 && in->eval_eps_exp <= 12) isog_eval_eps_exp = in->eval_eps_exp;
    isog_dform = in->dform;
    if (in->flood >= 0 && in->flood <= 2) isog_flood = in->flood;
}

const char* md_gto_iso_gpu_profile_variant(void) {
    static char buf[256];
    snprintf(buf, sizeof(buf), "%s | eval_eps 1e-%d dform %d flood %d", ISOG_REVISION, isog_eval_eps_exp, isog_dform, isog_flood);
    return buf;
}

// Waits for the GPU and appends a step: the time since the previous one. Returns its index (or ~0u).
static uint32_t isog_mark(md_gpu_stream_t s, const char* name, uint32_t items, uint32_t aux) {
    if (!isog_prof) return ~0u;
    md_gpu_stream_sync(s);
    const md_tick_t now = md_tick_now();
    const uint32_t i = isog_prof->num_steps;
    if (i < MD_GTO_ISO_GPU_PROFILE_CAP) {
        md_gto_iso_gpu_profile_step_t* st = &isog_prof->step[i];
        st->name = name;
        st->ms = md_tick_to_seconds(now - isog_prof_tick) * 1e3;
        st->items = items;
        st->aux = aux;
        isog_prof->num_steps = i + 1;
    }
    isog_prof_tick = now;
    return i < MD_GTO_ISO_GPU_PROFILE_CAP ? i : ~0u;
}

static void isog_mark_set(uint32_t i, uint32_t items, uint32_t aux) {
    if (!isog_prof || i >= isog_prof->num_steps) return;
    isog_prof->step[i].items = items;
    isog_prof->step[i].aux = aux;
}

// What happened since the previous mark is not timed (profile-only readbacks)
static void isog_mark_skip(md_gpu_stream_t s) {
    if (!isog_prof) return;
    md_gpu_stream_sync(s);
    isog_prof_tick = md_tick_now();
}

// BENCHMARK: kernel i is launched in this extraction (the profile lists the driver's statistics of it)
static void isog_used(int i) {
    if (isog_prof && i >= 0 && i < ISOG_K_COUNT) isog_prof->kernels_used |= 1ull << i;
}

static void isog_prof_begin(md_gpu_stream_t s) {
    if (!isog_prof) return;
    isog_prof->num_steps = 0;
    isog_prof->kernels_used = 0;
    isog_prof->brick_ao_sum = 0;
    isog_prof->brick_ao_max = 0;
    isog_prof->corner_cols = 0;
    isog_prof->corner_points = 0;
    isog_prof->cgroup_ao_sum = 0; isog_prof->cgroup_count = 0; isog_prof->cgroup_size = 0;
    isog_prof->vgroup_ao_sum = 0; isog_prof->vgroup_count = 0; isog_prof->vgroup_size = 0;
    isog_prof->brick_surf = 0; isog_prof->brick_in = 0; isog_prof->brick_out = 0;
    md_gpu_stream_sync(s);
    double best = 1e30;
    for (int i = 0; i < 8; ++i) {
        const md_tick_t t0 = md_tick_now();
        md_gpu_stream_sync(s);
        best = fmin(best, md_tick_to_seconds(md_tick_now() - t0) * 1e3);
    }
    isog_prof->wait_ms = best;
    isog_prof_tick = md_tick_now();
}

// Vertices still moving (vdone == 0) and the 64-vertex groups holding any of them (profile only)
static void isog_prof_active(md_gpu_stream_t s, md_gpu_addr_t vdone, uint32_t nv, uint32_t* out_active, uint32_t* out_groups) {
    *out_active = 0;
    *out_groups = 0;
    if (!isog_prof || nv == 0) return;
    uint32_t* d = (uint32_t*)malloc(sizeof(uint32_t) * nv);
    if (d && isog_read(s, d, vdone, sizeof(uint32_t) * nv)) {
        for (uint32_t g = 0; g < nv; g += 64) {
            uint32_t n = 0;
            for (uint32_t k = g; k < g + 64 && k < nv; ++k) n += d[k] == 0;
            *out_active += n;
            *out_groups += n > 0;
        }
    }
    free(d);
    isog_mark_skip(s);
}

size_t md_gto_iso_gpu_profile_print(const md_gto_iso_gpu_profile_t* p, char* buf, size_t cap) {
    if (!buf || cap == 0) return 0;
    buf[0] = '\0';
    if (!p) return 0;
    size_t len = 0;
#define ISOG_PR(...) do { if (len < cap) { const int n_ = snprintf(buf + len, cap - len, __VA_ARGS__); if (n_ > 0) len += (size_t)n_; } } while (0)
    enum { MAXN = 64 };
    const char* name[MAXN];
    uint32_t cnt[MAXN];
    double ms[MAXN];
    uint64_t items[MAXN], aux[MAXN];
    int nn = 0;
    double total = 0.0;
    for (uint32_t i = 0; i < p->num_steps; ++i) {
        const md_gto_iso_gpu_profile_step_t* st = &p->step[i];
        int k = 0;
        while (k < nn && strcmp(name[k], st->name) != 0) ++k;
        if (k == nn) {
            if (nn == MAXN) continue;
            name[k] = st->name; cnt[k] = 0; ms[k] = 0.0; items[k] = 0; aux[k] = 0;
            nn++;
        }
        cnt[k]++;
        ms[k] += st->ms;
        items[k] += st->items;
        aux[k] += st->aux;
        total += st->ms;
    }
    ISOG_PR("    %-26s %4s %9s %9s %11s %11s\n", "step", "n", "ms", "-waits", "items", "aux");
    for (int k = 0; k < nn; ++k) {
        const double w = ms[k] - cnt[k] * p->wait_ms;
        ISOG_PR("    %-26s %4u %9.3f %9.3f %11llu %11llu\n", name[k], cnt[k], ms[k], w > 0.0 ? w : 0.0, (unsigned long long)items[k], (unsigned long long)aux[k]);
    }
    ISOG_PR("    %-26s %4u %9.3f %9.3f   (one empty wait: %.3f ms)\n", "all", p->num_steps, total, fmax(0.0, total - p->num_steps * p->wait_ms), p->wait_ms);
    // work: the corner GEMM's multiply-adds (points x reaching AOs x function columns), the octree's
    for (int k = 0; k < nn; ++k) {
        const double w = fmax(1e-6, ms[k] - cnt[k] * p->wait_ms);
        if (strncmp(name[k], "corners (", 9) == 0 && items[k] > 0) {
            const double cols = (double)p->corner_cols;
            // the points evaluated per brick: 128 threads' worth, or those left after sharing (counted)
            const double pts = p->corner_points ? (double)p->corner_points / (double)items[k] : 128.0;
            // dform (no columns): A^2 / 2 per point for A AOs reaching the brick (the group's list is longer)
            const double fma = cols > 0.0 ? pts * (double)p->brick_ao_sum * cols : pts * 0.5 * (double)p->brick_ao_sq;
            ISOG_PR("    corners: %.1f AOs reach a brick (max %u of %u), %.1f points evaluated per brick, %u functions%s: %.3g multiply-adds, %.1f G/s\n",
                    (double)p->brick_ao_sum / (double)items[k], p->brick_ao_max, p->num_ao, pts, p->num_fn, cols > 0.0 ? "" : " (as D: A^2/2)",
                    fma, fma / (w * 1e6));
        }
        if (strncmp(name[k], "corners (", 9) == 0 && items[k] > 0 && p->brick_surf + p->brick_in + p->brick_out > 0) {
            ISOG_PR("    corner bricks: %u with surface cells, %u without (all corners inside %u, outside %u)\n",
                    p->brick_surf, p->brick_in + p->brick_out, p->brick_in, p->brick_out);
        }
        if (strncmp(name[k], "corners (", 9) == 0 && items[k] > 0 && p->cgroup_count > 0) {
            ISOG_PR("    corner groups: %.1f AOs reach a group of %u points (%u groups)\n",
                    (double)p->cgroup_ao_sum / (double)p->cgroup_count, p->cgroup_size, p->cgroup_count);
        }
        if (strcmp(name[k], "vertex eval") == 0 && p->vgroup_count > 0) {
            ISOG_PR("    vertex groups: %.1f AOs reach a group of %u vertices (%u groups)\n",
                    (double)p->vgroup_ao_sum / (double)p->vgroup_count, p->vgroup_size, p->vgroup_count);
        }
        if (strcmp(name[k], "octree level") == 0 && items[k] > 0 && (p->kernels_used & (1ull << ISOG_K_ODFORM))) {
            ISOG_PR("    octree (dform): %llu nodes, about 4 A^2 multiply-adds each for A AOs reaching a node\n", (unsigned long long)items[k]);
        } else if (strcmp(name[k], "octree level") == 0 && items[k] > 0) {
            const double fma = (double)items[k] * (double)p->num_ao * (double)p->fn_stride * 5.0;
            ISOG_PR("    octree: %llu nodes x %u AOs x %u function columns x 5: %.3g multiply-adds, %.1f G/s\n",
                    (unsigned long long)items[k], p->num_ao, p->fn_stride, fma, fma / (w * 1e6));
        }
    }
    // per pass: vertices still moving at each evaluation (and the groups holding them)
    ISOG_PR("    variant: %s\n", md_gto_iso_gpu_profile_variant());
    ISOG_PR("    eval passes, moving vertices (groups):");
    int col = 0;
    for (uint32_t i = 0; i < p->num_steps; ++i) {
        if (strcmp(p->step[i].name, "vertex eval") != 0) continue;
        if (col && col % 8 == 0) ISOG_PR("\n%42s", "");
        ISOG_PR(" %u (%u)", p->step[i].items, p->step[i].aux);
        col++;
    }
    ISOG_PR("\n");
    ISOG_PR("    octree levels, nodes in -> kept (ms):");
    for (uint32_t i = 0; i < p->num_steps; ++i) {
        if (strcmp(p->step[i].name, "octree level") != 0) continue;
        ISOG_PR(" %u->%u (%.2f)", p->step[i].items, p->step[i].aux, p->step[i].ms);
    }
    ISOG_PR("\n");
    // the main kernels as the driver compiled them (registers, shared memory: what limits occupancy)
    for (int i = 0; i < ISOG_K_COUNT; ++i) {
        md_gpu_kernel_info_t ki;
        if (!(p->kernels_used & (1ull << i)) || !isog_k[i] || !md_gpu_kernel_info(isog_k[i], &ki) || !ki.stats[0]) continue;
        ISOG_PR("    kernel %s (%u threads): %s\n", isog_kname[i] ? isog_kname[i] : "?", ki.group_size[0] * ki.group_size[1] * ki.group_size[2], ki.stats);
    }
#undef ISOG_PR
    return len < cap ? len : cap - 1;
}

// ---------------------------------------------------------------------------------------------
// The extraction as a job: the host steps through it, and where it needs a count from the device it
// copies the count to host-readable memory, records a sync and returns; md_gto_iso_gpu_job_poll takes
// the count once the sync has passed and issues the work after it.
// ---------------------------------------------------------------------------------------------

enum {
    ISOG_J_OCTREE,          // issue the octree (an attempt) and the bricks' slots, then read their counts
    ISOG_J_OCTREE_READ,     // the brick and slot counts: retry with more room, or on
    ISOG_J_BRICKS,          // the corners' buffers, the hash
    ISOG_J_PASS,            // a pass of the corners (one, or the flood's rounds)
    ISOG_J_FLOOD_READ,      // the bricks the next round adds
    ISOG_J_NETS,            // count the cells, then read the vertex count
    ISOG_J_NETS_READ,       // the vertex count: the vertices and the output, then read the final counts (and the mesh)
    ISOG_J_FINAL_READ,      // the final counts, the buckets: the mesh's ranges
    ISOG_J_MESH_READ,       // the rest of the mesh's host copy, when more quads came than the first copy held
    ISOG_J_DONE,
    ISOG_J_FAILED,
};

#define ISOG_JOB_ALLOCS 64
#define ISOG_JOB_FETCH  8
#define ISOG_JOB_MARKS  128

struct md_gto_iso_gpu_job_t {
    md_gpu_stream_t s;
    const md_gto_iso_gpu_set_t* set;
    md_gto_iso_gpu_desc_t desc;
    md_gto_iso_gpu_mesh_t out;
    int   stage;
    float progress;
    // device memory the job owns, freed at its end (recorded on the stream: the GPU may still use it)
    md_gpu_addr_t alloc[ISOG_JOB_ALLOCS];
    uint32_t num_alloc;
    // the readback in flight: copies into one host-readable buffer, then to their destinations
    bool          fetching, fetch_keep;
    md_gpu_mem_t  fetch_mem, kept;     // kept: a readback kept as the mesh's host copy
    size_t        mesh_off[3];         // in it: the positions, normals and indices
    uint32_t      idx_copied;          // indices it holds
    md_gpu_sync_t fetch_sync;
    uint32_t      num_fetch;
    void*         fetch_dst[ISOG_JOB_FETCH];
    size_t        fetch_off[ISOG_JOB_FETCH], fetch_bytes[ISOG_JOB_FETCH];
    // the extraction, across its waits
    md_tick_t t_start, ts;
    uint32_t F, root_mask;
    int L, Lt, l0, attempt;
    double origin[3];
    isog_field_t fields[MD_GTO_ISO_GPU_MAX_FIELDS];
    isog_args_t a;
    bool cdform, odform, flood;
    uint32_t cap, nb, ctr[12], cn[12], level_step[10];
    md_gpu_addr_t nodes[2], counts, bricks, oscan, bslot;
    md_gpu_uint4* bl;
    uint32_t* slot;
    uint32_t nslots, nkept, hcap;
    uint32_t pass, flood_rounds, flood_bricks, points_all, next;
    uint32_t nv, nq, cap_v, cap_q;
    uint32_t bk[128];
    // progress: markers recorded after chunks of the long launches, each the progress once it has passed
    md_gpu_sync_t mark_sync[ISOG_JOB_MARKS];
    float         mark_val[ISOG_JOB_MARKS];
    uint32_t      mark_head, mark_tail;
    uint32_t      chunks;          // the corner and vertex launches in this many chunks (1: whole)
    bool          blocking;        // md_gto_iso_gpu_extract's: nobody watches the progress, no chunks
};

static md_gpu_addr_t isog_jalloc(md_gto_iso_gpu_job_t* j, size_t bytes) {
    if (j->num_alloc >= ISOG_JOB_ALLOCS) return 0;
    const md_gpu_mem_t m = md_gpu_malloc(j->s, MD_GPU_MEM_DEVICE, bytes ? bytes : 16);
    if (m.gpu) j->alloc[j->num_alloc++] = m.gpu;
    return m.gpu;
}

static void isog_jfree(md_gto_iso_gpu_job_t* j, md_gpu_addr_t p) {
    for (uint32_t i = 0; i < j->num_alloc; ++i) {
        if (j->alloc[i] == p) {
            md_gpu_free(j->s, p);
            j->alloc[i] = j->alloc[--j->num_alloc];
            return;
        }
    }
}

// A copy of [src, src + bytes) to dst, made with the next isog_jfetch
static void isog_jfetch_add(md_gto_iso_gpu_job_t* j, void* dst, md_gpu_addr_t src, size_t bytes) {
    if (j->num_fetch >= ISOG_JOB_FETCH || !bytes) return;
    j->fetch_dst[j->num_fetch] = dst;
    j->fetch_off[j->num_fetch] = (size_t)src;      // the source until isog_jfetch
    j->fetch_bytes[j->num_fetch] = bytes;
    j->num_fetch++;
}

// Copies what was added to host-readable memory and records the sync the host waits for
static bool isog_jfetch(md_gto_iso_gpu_job_t* j) {
    size_t total = 0;
    for (uint32_t i = 0; i < j->num_fetch; ++i) total += (j->fetch_bytes[i] + 15) & ~(size_t)15;
    j->fetch_mem = md_gpu_malloc(j->s, MD_GPU_MEM_HOST_READ, total ? total : 16);
    if (!j->fetch_mem.gpu) { j->num_fetch = 0; return false; }
    size_t off = 0;
    for (uint32_t i = 0; i < j->num_fetch; ++i) {
        md_gpu_copy(j->s, j->fetch_mem.gpu + off, (md_gpu_addr_t)j->fetch_off[i], j->fetch_bytes[i]);
        j->fetch_off[i] = off;
        off += (j->fetch_bytes[i] + 15) & ~(size_t)15;
    }
    j->fetch_sync = md_gpu_stream_record(j->s);
    j->fetching = true;
    return true;
}

// The readback has passed: its bytes to their destinations
static void isog_jfetch_take(md_gto_iso_gpu_job_t* j) {
    for (uint32_t i = 0; i < j->num_fetch; ++i) {
        if (j->fetch_dst[i]) memcpy(j->fetch_dst[i], (const char*)j->fetch_mem.cpu + j->fetch_off[i], j->fetch_bytes[i]);
    }
    if (j->fetch_keep) j->kept = j->fetch_mem;
    else md_gpu_free(j->s, j->fetch_mem.gpu);
    j->fetch_mem = (md_gpu_mem_t){0};
    j->num_fetch = 0;
    j->fetching = false;
    j->fetch_keep = false;
}

static void isog_job_fail(md_gto_iso_gpu_job_t* j, const char* msg) {
    if (msg) MD_LOG_ERROR("%s", msg);
    md_gto_iso_gpu_mesh_free(j->s, &j->out);
    j->stage = ISOG_J_FAILED;
}

// The progress, never back
static void isog_jprog(md_gto_iso_gpu_job_t* j, float v) {
    if (v > j->progress) j->progress = v;
}

// A progress marker: the work issued so far is at 'value' once this sync has passed
static void isog_jmark(md_gto_iso_gpu_job_t* j, float value) {
    if (j->mark_tail - j->mark_head >= ISOG_JOB_MARKS) return;
    const uint32_t i = j->mark_tail++ % ISOG_JOB_MARKS;
    j->mark_sync[i] = md_gpu_stream_record(j->s);
    j->mark_val[i] = value;
}

// The markers passed: their progress
static void isog_jmarks_take(md_gto_iso_gpu_job_t* j) {
    while (j->mark_head != j->mark_tail) {
        const uint32_t i = j->mark_head % ISOG_JOB_MARKS;
        if (!md_gpu_sync_is_complete(j->mark_sync[i])) break;
        if (j->mark_val[i] > j->progress) j->progress = j->mark_val[i];
        j->mark_head++;
    }
}

// A launch of 'groups' groups in the job's chunks (direct grid), each followed by a marker, the progress
// going from p0 to p1 over them; whole when the job is not chunked. node_index adds group_base.
static void isog_jlaunch(md_gto_iso_gpu_job_t* j, md_gpu_kernel_t k, isog_args_t args, uint32_t groups, float p0, float p1) {
    const uint32_t n = j->chunks > 1 ? j->chunks : 1;
    uint32_t per = (groups + n - 1) / n;
    if (per > 65535u) per = 65535u;     // one row of groups a chunk: the grid's y adds no groups past it
    if (n == 1 || per == 0) {
        MD_GPU_LAUNCH(j->s, k, isog_groups(groups), args);
        return;
    }
    const uint32_t m = (groups + per - 1) / per;
    for (uint32_t c = 0; c < m; ++c) {
        args.group_base = c * per;
        const uint32_t g = groups - c * per < per ? groups - c * per : per;
        MD_GPU_LAUNCH(j->s, k, md_gpu_grid(g, 1, 1), args);
        isog_jmark(j, p0 + (p1 - p0) * (float)(c + 1) / (float)m);
    }
}

// The same for a launch whose groups the device counts (grid from count_in by iso_grid): 'expect' groups
// expected, the last chunk taking whatever is left
static void isog_jlaunch_indirect(md_gto_iso_gpu_job_t* j, md_gpu_kernel_t k, isog_args_t args, uint32_t expect, float p0, float p1) {
    md_gpu_stream_t s = j->s;
    const uint32_t n = j->chunks > 1 ? j->chunks : 1;
    uint32_t per = (expect + n - 1) / n;
    if (per > 65535u) per = 65535u;
    if (n == 1 || per == 0) {
        MD_GPU_LAUNCH(s, isog_k[ISOG_K_GRID], md_gpu_grid(1, 1, 1), args);
        md_gpu_launch_indirect(s, k, args.grid_out, &args, sizeof(args));
        return;
    }
    const uint32_t m = (expect + per - 1) / per;
    for (uint32_t c = 0; c < m; ++c) {
        args.group_base = c * per;
        args.grid_max = c + 1 < m ? per : 0xFFFFFFFFu;     // the last: the rest, however many
        MD_GPU_LAUNCH(s, isog_k[ISOG_K_GRID], md_gpu_grid(1, 1, 1), args);
        md_gpu_launch_indirect(s, k, args.grid_out, &args, sizeof(args));
        isog_jmark(j, p0 + (p1 - p0) * (float)(c + 1) / (float)m);
    }
}

// Runs the job from its stage until it waits for a readback, is done, or fails
static void isog_job_step(md_gto_iso_gpu_job_t* j) {
    md_gpu_stream_t s = j->s;
    const md_gto_iso_gpu_set_t* set = j->set;
    const md_gto_iso_gpu_desc_t* desc = &j->desc;
    md_gto_iso_gpu_mesh_t* out = &j->out;
    const uint32_t F = j->F;
    const int L = j->L;
    const double h = desc->cell_size;
    isog_args_t* ap0 = &j->a;
#define a (*ap0)
    while (!j->fetching && j->stage != ISOG_J_DONE && j->stage != ISOG_J_FAILED) {
        switch (j->stage) {
        case ISOG_J_OCTREE: {
            // ---- octree, one launch per level sized by the previous level's count; retried with more room
            const uint32_t cap = j->cap;
            j->nodes[0] = isog_jalloc(j, sizeof(md_gpu_uint4) * cap);
            j->nodes[1] = isog_jalloc(j, sizeof(md_gpu_uint4) * cap);
            j->counts = isog_jalloc(j, sizeof(uint32_t) * (L + 2));
            j->oscan = isog_jalloc(j, sizeof(uint32_t) * ((size_t)cap + 1));     // a level's counts, then offsets
            j->bslot = isog_jalloc(j, sizeof(uint32_t) * ((size_t)cap + 1));     // the bricks' first slots
            if (!j->nodes[0] || !j->nodes[1] || !j->counts || !j->oscan || !j->bslot) { isog_job_fail(j, "md_gto_iso_gpu_extract: octree overflow or out of device memory"); break; }
            md_gpu_memset(s, j->counts, 0, sizeof(uint32_t) * (L + 2));
            md_gpu_memset(s, a.counters, 0, sizeof(uint32_t) * 12);
            // The bricks' level is not tested: every child of a node of the last tested level is a brick, with
            // its parent's sheet bits, and the corners decide (for less than the bound costs). A superset of
            // the bricks, so no cell is missed.
            const int Lt = L - 1;
            // the first tested level l0, its nodes every node of the levels above (siblings together, in
            // Morton order), with the root's sheet bits
            const int l0 = Lt > 0 ? (ISOG_OCTREE_START < Lt ? ISOG_OCTREE_START : Lt) : 0;
            j->Lt = Lt;
            j->l0 = l0;
            {
                md_gpu_uint4 lv[512];
                uint32_t nl = 1;
                const uint32_t root_mask = j->root_mask;
                lv[0] = (md_gpu_uint4){0, 0, 0, root_mask};
                for (int l = 0; l < l0; ++l) {
                    md_gpu_uint4 nx[512];
                    for (uint32_t i = 0; i < nl; ++i)
                        for (uint32_t c = 0; c < 8; ++c)
                            nx[8 * i + c] = (md_gpu_uint4){lv[i].x * 2 + (c & 1), lv[i].y * 2 + ((c >> 1) & 1), lv[i].z * 2 + (c >> 2), root_mask};
                    nl *= 8;
                    memcpy(lv, nx, sizeof(md_gpu_uint4) * nl);
                }
                md_gpu_upload(s, j->nodes[l0 & 1], lv, sizeof(md_gpu_uint4) * nl);
                md_gpu_upload(s, j->counts + sizeof(uint32_t) * l0, &nl, sizeof(nl));
            }
            isog_mark(s, "octree setup", 0, 0);
            bool ok = true;
            for (int l = l0; l <= Lt; ++l) {
                isog_args_t o = a;
                o.level = (uint32_t)l;
                o.nodes_in  = j->nodes[l & 1];
                o.nodes_out = j->nodes[(l + 1) & 1];
                o.count_in  = j->counts + sizeof(uint32_t) * l;
                o.count_out = j->counts + sizeof(uint32_t) * (l + 1);
                o.cap = cap;
                o.grid_mul = 1;
                o.grid_div = 0;
                o.mode = 0u;
                o.scan = j->oscan;
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_GRID], md_gpu_grid(1, 1, 1), o);
                // the D form's bound for densities; a thread per shell for at most 4 functions; else the
                // functions' ranges in groupshared memory
                int ko = j->odform ? ISOG_K_ODFORM : set->num_fn <= 4 ? ISOG_K_OCTREE_SMALL : ISOG_K_OCTREE;
                if (ko == ISOG_K_OCTREE_SMALL && !isog_k[ko]) ko = ISOG_K_OCTREE;
                if (!isog_have(ko)) { ok = false; break; }
                md_gpu_launch_indirect(s, isog_k[ko], a.grid_out, &o, sizeof(o));
                isog_used(ko);
                {   // in order: the offsets of what the nodes emit (their prefix sum), then the emission
                    isog_args_t e = o;
                    e.count = 0xFFFFFFFFu;     // as many as count_in says; the sum to count_out
                    MD_GPU_LAUNCH(s, isog_k[ISOG_K_SCAN], md_gpu_grid(1, 1, 1), e);
                    e.grid_div = 64;
                    MD_GPU_LAUNCH(s, isog_k[ISOG_K_GRID], md_gpu_grid(1, 1, 1), e);
                    md_gpu_launch_indirect(s, isog_k[ISOG_K_OEMIT], a.grid_out, &e, sizeof(e));
                }
                j->level_step[l] = isog_mark(s, "octree level", 0, 0);
            }
            if (!ok) { isog_job_fail(j, "md_gto_iso_gpu_extract: octree overflow or out of device memory"); break; }
            if (isog_prof) {    // nodes per level: in, and kept (children appended, 8 per kept node above the bricks)
                uint32_t cn[12] = {0};
                isog_read(s, cn, j->counts, sizeof(uint32_t) * (L + 2));
                for (int l = l0; l <= Lt; ++l) isog_mark_set(j->level_step[l], cn[l], l < L ? cn[l + 1] / 8 : cn[l + 1]);
                isog_mark_skip(s);
            }
            j->bricks = j->nodes[(Lt + 1) & 1];
            {   // ---- the bricks' corner slots: one per (brick, field with a sheet bit), their prefix sum
                isog_args_t b = a;
                b.bricks = j->bricks;
                b.count_in = j->counts + sizeof(uint32_t) * (Lt + 1);
                b.count_out = a.counters + sizeof(uint32_t) * 4;      // the slots
                b.cap = cap;
                b.scan = j->bslot;
                b.grid_mul = 1;
                b.grid_div = 64;
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_GRID], md_gpu_grid(1, 1, 1), b);
                md_gpu_launch_indirect(s, isog_k[ISOG_K_BSLOTS], a.grid_out, &b, sizeof(b));
                b.count = 0xFFFFFFFFu;
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_SCAN], md_gpu_grid(1, 1, 1), b);
                isog_mark(s, "brick slots", 0, 0);
            }
            isog_jfetch_add(j, &j->nb, j->counts + sizeof(uint32_t) * (Lt + 1), sizeof(uint32_t));
            isog_jfetch_add(j, j->ctr, a.counters, sizeof(j->ctr));
            if (!isog_jfetch(j)) { isog_job_fail(j, "md_gto_iso_gpu_extract: out of device memory"); break; }
            j->stage = ISOG_J_OCTREE_READ;
            isog_jprog(j, 0.05f);
        } break;

        case ISOG_J_OCTREE_READ: {
            isog_mark(s, "readback: brick count", j->nb, 0);
            if (j->ctr[3] || j->nb > j->cap) {
                if (++j->attempt < 4) {     // more room, and again
                    isog_jfree(j, j->nodes[0]);
                    isog_jfree(j, j->nodes[1]);
                    isog_jfree(j, j->counts);
                    isog_jfree(j, j->oscan);
                    isog_jfree(j, j->bslot);
                    j->cap *= 4;
                    j->nb = 0;
                    j->stage = ISOG_J_OCTREE;
                    break;
                }
                isog_job_fail(j, "md_gto_iso_gpu_extract: octree overflow or out of device memory");
                break;
            }
            isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_OCTREE, &j->ts);
            // A long octree means a long extraction: its corner and vertex launches then go in chunks, each
            // followed by a progress marker (and a point where other work, the rendering, can get the GPU).
            // The octree is a tenth to a quarter of an extraction: a chunk every 2 ms of it, about 10-20 ms of
            // work each. A quick one (under 4 ms) is not split: the markers would cost more than they tell.
            {
                const double oms = md_tick_to_seconds(md_tick_now() - j->t_start) * 1e3;
                uint32_t n = (oms < 4.0 || j->blocking || isog_prof) ? 1u : (uint32_t)(oms / 2.0);
                if (n > 32) n = 32;
                if (isog_chunks > 0) n = (uint32_t)isog_chunks;
                j->chunks = n;
            }
            // the bricks in Morton order and their slots, all on the device; the counts sized what follows
            j->nslots = j->ctr[4];
            j->nkept = j->ctr[8];     // the bricks may keep their sibling slots: rejected ones have no sheet bits
            out->num_bricks = j->nkept;
            out->num_corner_slots = j->nslots;
            j->stage = ISOG_J_BRICKS;
            isog_jprog(j, 0.15f);
        } break;

        case ISOG_J_BRICKS: {
            const uint32_t nb = j->nb, nslots = j->nslots;
            const double* origin = j->origin;
            if (isog_prof || isog_debug) {    // the bricks and their slots on the host: the profile's and the tests' only
                j->bl = (md_gpu_uint4*)malloc(sizeof(md_gpu_uint4) * (nb + 1));
                j->slot = (uint32_t*)malloc(sizeof(uint32_t) * (nb + 1));
                if (j->bl && j->slot && nb) {
                    isog_read(s, j->bl, j->bricks, sizeof(md_gpu_uint4) * nb);
                    isog_read(s, j->slot, j->bslot, sizeof(uint32_t) * nb);
                }
                isog_mark_skip(s);
            }
            const md_gpu_uint4* bl = j->bl;
            if (isog_prof && bl) {    // BENCHMARK: the AOs that reach each brick (what the corner kernels list), not timed
                uint64_t sum = 0, sq = 0;
                uint32_t mx = 0;
                const float he = (float)h;
                for (uint32_t b = 0; b < nb; ++b) {
                    if (bl[b].w == 0) continue;
                    const float blo[3] = {(float)origin[0] + (float)(bl[b].x * 4u) * he, (float)origin[1] + (float)(bl[b].y * 4u) * he, (float)origin[2] + (float)(bl[b].z * 4u) * he};
                    const float bhi[3] = {blo[0] + 4.0f * he, blo[1] + 4.0f * he, blo[2] + 4.0f * he};
                    const uint32_t n = isog_prof_box_aos(set, blo, bhi);
                    sum += n;
                    sq += (uint64_t)n * n;
                    mx = n > mx ? n : mx;
                }
                isog_prof->brick_ao_sum = sum;
                isog_prof->brick_ao_sq = sq;
                isog_prof->brick_ao_max = mx;
                isog_prof->num_ao = set->num_ao;
                isog_prof->num_fn = set->num_fn;
                isog_prof->fn_stride = set->fn_stride;
                isog_prof->orbitals_only = set->orbitals_only;
                isog_mark_skip(s);
            }
            a.bricks = j->bricks;
            a.count = nb;
            a.brick_slot = j->bslot;
            uint32_t hcap = 64;
            while (hcap < 2 * nb) hcap <<= 1;
            j->hcap = hcap;
            a.hash_mask = hcap - 1;
            a.hash_keys = isog_jalloc(j, sizeof(uint32_t) * hcap);
            a.hash_vals = isog_jalloc(j, sizeof(uint32_t) * hcap);
            a.corners   = isog_jalloc(j, sizeof(float) * 125 * ((size_t)nslots + 1));
            a.cell_vert = isog_jalloc(j, sizeof(uint32_t) * 128 * ((size_t)nslots + 1));
            a.plist = isog_jalloc(j, sizeof(uint32_t) * 125 * ((size_t)nb + 1));
            // per brick: counts, then offsets (ordered lists). Allocated either way: llvmpipe may load it in a
            // branch not taken
            a.scan  = isog_jalloc(j, sizeof(uint32_t) * ((size_t)nb + 1));
            a.bstate = isog_jalloc(j, sizeof(uint32_t) * ((size_t)nb + 1));     // flood: per brick, its round
            if (!a.scan || !a.bstate || !a.brick_slot || !a.hash_keys || !a.hash_vals || !a.corners || !a.cell_vert || !a.plist) {
                isog_job_fail(j, "md_gto_iso_gpu_extract: out of device memory or buffer overflow");
                break;
            }
            isog_jprog(j, 0.2f);
            if (nb == 0) {
                isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_CORNERS, &j->ts);
                j->stage = ISOG_J_NETS;
                break;
            }
            md_gpu_memset(s, a.hash_keys, 0, sizeof(uint32_t) * hcap);
            MD_GPU_LAUNCH(s, isog_k[ISOG_K_HASH], md_gpu_grid_for(isog_k[ISOG_K_HASH], nb, 1, 1), a);
            isog_mark(s, "hash", nb, hcap);
            // the points to evaluate, listed brick by brick; then 256 of them per group. APPROXIMATE: seed and
            // flood (GTO_ISO_FLOOD): pass 0 the bricks' corners, then the seeds (and the bricks holding an atom),
            // then round by round the bricks across faces with a sign change; else one pass over every brick
            j->pass = 0;
            j->flood_rounds = j->flood_bricks = j->points_all = 0;
            j->stage = ISOG_J_PASS;
        } break;

        case ISOG_J_PASS: {
            const uint32_t nb = j->nb, pass = j->pass, nfn = set->num_fn;
            const bool flood = j->flood, cdform = j->cdform;
            const uint32_t nkept = j->nkept, nslots = j->nslots;
            const double* origin = j->origin;
            const md_gpu_uint4* bl = j->bl;
            md_gpu_memset(s, a.counters + sizeof(uint32_t) * 5, 0, sizeof(uint32_t));
            isog_args_t ap = a;
            // shared corners (mode 1): a brick's high-face corners from the neighbour that owns them
            ap.mode = 1u;
            if (flood) { ap.mode |= pass == 0 ? 16u : 8u; ap.level = pass; }
            // in the order of the bricks: the counts, their prefix sum, then the points at the offsets
            ap.mode |= 2u;
            MD_GPU_LAUNCH(s, isog_k[ISOG_K_POINTS], isog_groups(nb), ap);
            MD_GPU_LAUNCH(s, isog_k[ISOG_K_SCAN], md_gpu_grid(1, 1, 1), ap);
            isog_used(ISOG_K_SCAN);
            ap.mode = (ap.mode & ~2u) | 4u;
            MD_GPU_LAUNCH(s, isog_k[ISOG_K_POINTS], isog_groups(nb), ap);
            isog_mark(s, "corners: point list", nb, 0);
            // one function: four points per thread (iso_corners_list_f1v4); densities as matrices: iso_dcorners,
            // 128 points per group; canonical sets of more than 64 functions: the wide GEMM, 64 points per group;
            // else the tiled GEMM, 128 points per group (iso_gcorners_f*)
            const bool one = nfn == 1 && F == 1;
            const int kw = nfn <= 128 ? ISOG_K_GCORN_W128 : nfn <= 192 ? ISOG_K_GCORN_W192 : ISOG_K_GCORN_Y256;
            const bool wide = !one && !cdform && a.canon && nfn > 64u && isog_k[kw];
            const bool gemm = !one;
            const int kc = one ? ISOG_K_CLIST_F1V4 : cdform ? ISOG_K_DCORN : wide ? kw
                         : nfn <= 16 ? ISOG_K_GCORN_F16 : nfn <= 32 ? ISOG_K_GCORN_F32 : ISOG_K_GCORN_F64;
            if (!isog_have(kc)) { isog_job_fail(j, "md_gto_iso_gpu_extract: a corner kernel is missing"); break; }
            isog_args_t ac = a;
            ac.count_in = a.counters + sizeof(uint32_t) * 5;
            ac.cap = 0;
            ac.grid_mul = 1;
            ac.grid_div = wide ? 64 : gemm ? 128 : 256;
            {   // in chunks when the job is long: the points expected (about 70 a brick evaluated in full, as the
                // profile counts them; a brick's corners in the flood's first pass), the progress spread over them
                const uint32_t nbp = !flood ? nkept : pass == 0 ? nb : pass == 1 ? j->next : j->next;
                const uint64_t pts = (uint64_t)nbp * (flood && pass == 0 ? 8u : 70u);
                const uint32_t expect = (uint32_t)((pts + ac.grid_div - 1) / ac.grid_div);
                const float p0 = j->progress < 0.15f ? 0.15f : j->progress;
                const float p1 = !flood ? 0.6f : pass == 0 ? 0.2f : pass == 1 ? 0.55f : p0 + 0.5f * (0.6f - p0);
                isog_jlaunch_indirect(j, isog_k[kc], ac, (flood && pass == 0) ? 0u : expect, p0, p1);
            }
            isog_used(kc);
            uint32_t fc = 0;
            const char* kn = "";
            switch (kc) {
                case ISOG_K_CLIST_F1V4: fc = 1;  kn = "corners (list f1 x4)"; break;
                case ISOG_K_GCORN_F16:  fc = 16; kn = "corners (gemm f16)"; break;
                case ISOG_K_GCORN_F32:  fc = 32; kn = "corners (gemm f32)"; break;
                case ISOG_K_DCORN:     fc = 0;  kn = "corners (dform)"; break;
                case ISOG_K_GCORN_W128: fc = 128; kn = "corners (gemm w128)"; break;
                case ISOG_K_GCORN_W192: fc = 192; kn = "corners (gemm w192)"; break;
                case ISOG_K_GCORN_Y256: fc = 256; kn = "corners (gemm y256)"; break;
                default:                fc = 64; kn = "corners (gemm f64)"; break;
            }
            if (isog_prof) isog_prof->corner_cols = fc ? fc * ((nfn + fc - 1) / fc) : 0;
            const uint32_t ck = isog_mark(s, kn, nkept, nslots);
            if (isog_prof && flood) {   // aux: the lattice points evaluated, summed over the passes
                uint32_t c8[8];
                isog_read(s, c8, a.counters, sizeof(c8));
                j->points_all += c8[5];
                isog_mark_set(ck, nkept, c8[5]);
                isog_prof->corner_points = j->points_all;
                isog_mark_skip(s);
            } else if (isog_prof) {    // aux: the lattice points evaluated
                uint32_t c8[8];
                isog_read(s, c8, a.counters, sizeof(c8));
                isog_mark_set(ck, nkept, c8[5]);
                isog_prof->corner_points = c8[5];
                // the boxes of the groups the corner kernel takes (consecutive points of the list)
                const uint32_t np = c8[5];
                uint32_t* pl = (uint32_t*)malloc(sizeof(uint32_t) * (np + 1));
                float* px = (float*)malloc(sizeof(float) * 3 * (np + 1));
                if (pl && px && np > 0 && isog_read(s, pl, a.plist, sizeof(uint32_t) * np)) {
                    for (uint32_t i = 0; i < np; ++i) {
                        const uint32_t b = pl[i] >> 7, k = (pl[i] & 127u) < 124u ? (pl[i] & 127u) : 124u;
                        const uint32_t g[3] = {bl[b].x * 4u + k % 5u, bl[b].y * 4u + (k / 5u) % 5u, bl[b].z * 4u + k / 25u};
                        for (int ax = 0; ax < 3; ++ax) px[3 * i + ax] = (float)origin[ax] + (float)g[ax] * (float)h;
                    }
                    isog_prof->cgroup_size = wide ? 64u : gemm ? 128u : 256u;     // f1v4: 64 threads x 4 points
                    isog_prof_groups(set, px, 3, np, isog_prof->cgroup_size, &isog_prof->cgroup_ao_sum, &isog_prof->cgroup_count);
                }
                free(pl);
                free(px);
                isog_mark_skip(s);
            }
            {   // the corners a brick shares with the neighbours that evaluated them
                isog_args_t af = a;
                af.mode = flood ? (pass == 0 ? 16u : 8u) : 0u;
                af.level = pass;
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_CORNERS_FILL], isog_groups(nb), af);
                isog_mark(s, "corners: fill shared", nb, 0);
            }
            if (!flood) {
                isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_CORNERS, &j->ts);
                j->stage = ISOG_J_NETS;
                break;
            }
            // the next round: the seeds (pass 0), or the bricks across this round's faces with a sign change
            md_gpu_memset(s, a.counters + sizeof(uint32_t) * 6, 0, sizeof(uint32_t));
            isog_args_t ar = a;
            ar.level = pass;
            if (pass == 0) {
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_FSEED], md_gpu_grid((nb + 63) / 64, 1, 1), ar);
                isog_used(ISOG_K_FSEED);
                if (isog_flood_atoms) {
                    MD_GPU_LAUNCH(s, isog_k[ISOG_K_FATOMS], md_gpu_grid((set->num_shells + 63) / 64, 1, 1), ar);
                    isog_used(ISOG_K_FATOMS);
                }
            } else {
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_FFACE], isog_groups(nb), ar);
                isog_used(ISOG_K_FFACE);
            }
            j->next = 0;
            isog_jfetch_add(j, &j->next, a.counters + sizeof(uint32_t) * 6, sizeof(uint32_t));
            if (!isog_jfetch(j)) { isog_job_fail(j, "md_gto_iso_gpu_extract: out of device memory"); break; }
            j->stage = ISOG_J_FLOOD_READ;
        } break;

        case ISOG_J_FLOOD_READ: {
            const uint32_t nb = j->nb;
            isog_mark(s, j->pass == 0 ? "flood: seeds" : "flood: faces", nb, j->next);
            j->flood_bricks += j->next;
            if (j->next != 0 && j->pass < 64) {
                // the pass done: the end of its share (the first its brick corners, the second the most, the
                // others less and less)
                isog_jprog(j, j->pass == 0 ? 0.2f : j->pass == 1 ? 0.55f : j->progress + 0.5f * (0.6f - j->progress));
                j->flood_rounds++;
                j->pass++;
                j->stage = ISOG_J_PASS;
                break;
            }
            // the bricks no round reached: no cells (and the host's list says so, for the profile)
            MD_GPU_LAUNCH(s, isog_k[ISOG_K_FFINISH], md_gpu_grid((nb + 63) / 64, 1, 1), a);
            isog_used(ISOG_K_FFINISH);
            if (isog_prof || isog_debug) isog_read(s, j->bl, a.bricks, sizeof(md_gpu_uint4) * nb);   // the profile's and the tests' only
            isog_mark(s, "flood: finish", j->flood_bricks, j->flood_rounds);
            isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_CORNERS, &j->ts);
            j->stage = ISOG_J_NETS;
        } break;

        case ISOG_J_NETS: {
            // ---- Surface Nets: count, then size the vertex buffers and emit
            const uint32_t nb = j->nb;
            j->nv = 0;
            j->stage = ISOG_J_NETS_READ;
            if (nb == 0) break;
            a.mode = 0u;     // the bricks' counts into scan, for vertices in the order of the bricks
            md_gpu_memset(s, a.counters, 0, sizeof(uint32_t) * 7);     // [7] (dform overflow) kept
            MD_GPU_LAUNCH(s, isog_k[ISOG_K_CELLS], isog_groups(nb), a);
            MD_GPU_LAUNCH(s, isog_k[ISOG_K_SCAN], md_gpu_grid(1, 1, 1), a);
            isog_used(ISOG_K_SCAN);
            isog_mark(s, "cells: count", nb, 0);
            isog_jfetch_add(j, &j->nv, a.counters, sizeof(uint32_t));
            if (!isog_jfetch(j)) { isog_job_fail(j, "md_gto_iso_gpu_extract: out of device memory"); break; }
        } break;

        case ISOG_J_NETS_READ: {
            const uint32_t nb = j->nb, nslots = j->nslots, nv = j->nv;
            const md_gpu_uint4* bl = j->bl;
            const uint32_t* slot = j->slot;
            const isog_field_t* fields = j->fields;
            isog_jprog(j, 0.65f);
            if (nb > 0) isog_mark(s, "readback: vertex count", nv, 0);
            if (nb > 0 && isog_prof) {     // BENCHMARK: the bricks given corners, by their corners (not timed)
                uint32_t* sc = (uint32_t*)malloc(sizeof(uint32_t) * (nb + 1));
                float* cv = (float*)malloc(sizeof(float) * 125 * ((size_t)nslots + 1));
                if (sc && cv && isog_read(s, sc, a.scan, sizeof(uint32_t) * (nb + 1)) &&
                    isog_read(s, cv, a.corners, sizeof(float) * 125 * (size_t)nslots)) {
                    for (uint32_t b = 0; b < nb; ++b) {
                        const uint32_t m = bl[b].w;
                        if (m == 0) continue;
                        if (sc[b + 1] > sc[b]) { isog_prof->brick_surf++; continue; }
                        bool in = false;
                        for (uint32_t f = 0; f < F && !in; ++f) {
                            const uint32_t fb = (m >> (2 * f)) & 3u;
                            if (fb == 0) continue;
                            const float* c = cv + (size_t)(slot[b] + popcnt32((m | (m >> 1)) & 0x55555555u & ((1u << (2 * f)) - 1u))) * 125;
                            for (uint32_t sh = 0; sh < 2 && !in; ++sh) {
                                if (!(fb & (1u << sh))) continue;
                                const float sg = sh == 0 ? 1.0f : -1.0f;
                                bool all = true;
                                for (int k = 0; k < 125 && all; ++k) all = sg * c[k] > fields[f].iso;
                                in = all;
                            }
                        }
                        if (in) isog_prof->brick_in++; else isog_prof->brick_out++;
                    }
                }
                free(sc);
                free(cv);
                isog_mark_skip(s);
            }
            if (isog_debug) {
                isog_debug_t* dbg = isog_debug;
                dbg->num_bricks = nb;
                dbg->num_slots = nslots;
                dbg->bricks  = (uint32_t*)malloc(sizeof(md_gpu_uint4) * (nb + 1));
                dbg->slots   = (uint32_t*)malloc(sizeof(uint32_t) * (nb + 1));
                dbg->corners = (float*)malloc(sizeof(float) * 125 * ((size_t)nslots + 1));
                memcpy(dbg->bricks, bl, sizeof(md_gpu_uint4) * nb);
                memcpy(dbg->slots, slot, sizeof(uint32_t) * nb);
                isog_read(s, dbg->corners, a.corners, sizeof(float) * 125 * (size_t)nslots);
                memcpy(dbg->origin, j->origin, sizeof(j->origin));
                dbg->brick_level = L;
                isog_mark_skip(s);
            }
            free(j->bl);
            free(j->slot);
            j->bl = NULL;
            j->slot = NULL;

            const uint32_t cap_v = nv + 1, cap_q = 3 * nv + 1;
            j->cap_v = cap_v;
            j->cap_q = cap_q;
            a.vpos     = isog_jalloc(j, sizeof(md_gpu_float4) * cap_v);
            a.vpos_out = isog_jalloc(j, sizeof(md_gpu_float4) * cap_v);
            a.vnrm     = isog_jalloc(j, sizeof(md_gpu_float4) * cap_v);
            a.vval     = isog_jalloc(j, sizeof(md_gpu_float4) * cap_v);
            a.vinfo    = isog_jalloc(j, sizeof(md_gpu_uint4) * cap_v);
            a.quads    = isog_jalloc(j, sizeof(md_gpu_uint4) * cap_q);
            a.remap    = isog_jalloc(j, sizeof(uint32_t) * cap_v);
            a.bucket   = isog_jalloc(j, sizeof(uint32_t) * 192);
            a.vdone    = isog_jalloc(j, sizeof(uint32_t) * cap_v);
            if (!a.vpos || !a.vpos_out || !a.vnrm || !a.vval || !a.vinfo || !a.quads || !a.remap || !a.bucket || !a.vdone) {
                isog_job_fail(j, "md_gto_iso_gpu_extract: out of device memory or buffer overflow");
                break;
            }
            if (nv > 0) {
                md_gpu_memset(s, a.cell_vert, 0xFF, sizeof(uint32_t) * 128 * (size_t)nslots);
                md_gpu_memset(s, a.counters, 0, sizeof(uint32_t) * 7);     // [7] (dform overflow) kept
                a.mode = 1u;
                a.cap = cap_v;
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_CELLS], isog_groups(nb), a);
                isog_mark(s, "cells: emit", nb, nv);
                a.cap = cap_q;
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_QUADS], isog_groups(nb), a);
                isog_mark(s, "quads", nb, 0);
            }
            isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_NETS, &j->ts);

            // ---- vertices: relaxed on the lattice's own (trilinear) surface, with its normals from iso_cells,
            // then projected by Newton steps: no field evaluation before the projection
            if (nv > 0) {
                isog_args_t v = a;
                v.count = nv;
                const md_gpu_grid_t g_vert = md_gpu_grid_for(isog_k[ISOG_K_NEWTON], nv, 1, 1);
                // one orbital: the projection in one launch, a group of 128 vertices through every iteration
                // (iso_project_list_f1); other sets of orbitals: a thread per vertex (iso_eval_orbital); densities
                // (and sets with any): the tiled GEMM, 32 vertices (4 channels each) per group, the wide one for
                // canonical sets of more than 64 functions, 16 per group; or the D form
                const isog_vkernel_t vk = isog_vertex_kernel(set, a.dform != 0u, true);
                const bool one = vk.k == ISOG_K_PROJ1;
                if (!isog_have(vk.k)) { isog_job_fail(j, "md_gto_iso_gpu_extract: an evaluation kernel is missing"); break; }
                md_gpu_kernel_t k_eval = isog_k[vk.k];
                isog_used(vk.k);
                // the evaluation's groups (0: a thread per vertex)
                const uint32_t n_eval_groups = vk.per_group ? (nv + vk.per_group - 1) / vk.per_group : 0u;
                const md_gpu_grid_t g_eval = n_eval_groups ? isog_groups(n_eval_groups) : g_vert;
                if (isog_prof) {    // BENCHMARK: the boxes of the vertex groups (consecutive vertices), not timed
                    const uint32_t vg = vk.per_group;
                    md_gpu_float4* vp = vg ? (md_gpu_float4*)malloc(sizeof(md_gpu_float4) * nv) : NULL;
                    if (vp && isog_read(s, vp, a.vpos, sizeof(md_gpu_float4) * nv)) {
                        isog_prof->vgroup_size = vg;
                        isog_prof_groups(set, &vp[0].x, 4, nv, vg, &isog_prof->vgroup_ao_sum, &isog_prof->vgroup_count);
                    }
                    free(vp);
                    isog_mark_skip(s);
                }
                // relaxation reads vpos and writes vpos_out: the two swap roles after every pass (no copy)
                for (uint32_t r = 0; r < desc->relax_iterations; ++r) {
                    MD_GPU_LAUNCH(s, isog_k[ISOG_K_RELAX], g_vert, v);
                    isog_mark(s, "relax", nv, 0);
                    const md_gpu_addr_t t = v.vpos;
                    v.vpos = v.vpos_out;
                    v.vpos_out = t;
                }
                isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_RELAX, &j->ts);
                const uint32_t iters = desc->project_iterations;
                if (one) {
                    v.level = iters;
                    isog_jlaunch(j, k_eval, v, n_eval_groups, 0.65f, 0.92f);
                    isog_mark(s, "vertex project (fused)", nv, 0);
                } else {
                    // an evaluation and a step per iteration; converged vertices are flagged and skipped (they no
                    // longer move, so their values stay valid for the final normals)
                    md_gpu_memset(s, v.vdone, 0, sizeof(uint32_t) * nv);
                    for (uint32_t it = 0; it <= iters; ++it) {
                        uint32_t act, grp;
                        isog_prof_active(s, v.vdone, nv, &act, &grp);
                        if (it == 0 && n_eval_groups) isog_jlaunch(j, k_eval, v, n_eval_groups, 0.65f, 0.9f);
                        else MD_GPU_LAUNCH(s, k_eval, g_eval, v);
                        isog_mark(s, "vertex eval", act, grp);
                        v.mode = it == iters ? 1u : 0u;     // 1: the normals; 0: a step (one predicted far below the tolerance ends the vertex)
                        MD_GPU_LAUNCH(s, isog_k[ISOG_K_NEWTON], g_vert, v);
                        isog_mark(s, v.mode == 1u ? "vertex normals" : "vertex newton", act, 0);
                    }
                }
                isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_PROJECT, &j->ts);
                a.vpos = v.vpos;            // where the vertices ended up
                a.vpos_out = v.vpos_out;
            }
            // ---- regroup per (field, sheet) into the output buffers: sized by the vertices, the quads (on the
            // device: no count is read for them) at most cap_q
            md_gpu_mem_t pos = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(md_gpu_float4) * cap_v);
            md_gpu_mem_t nrm = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(md_gpu_float4) * cap_v);
            md_gpu_mem_t idx = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(uint32_t) * 6 * (size_t)cap_q);
            out->position = pos.gpu; out->normal = nrm.gpu; out->index = idx.gpu;
            if (!pos.gpu || !nrm.gpu || !idx.gpu) { isog_job_fail(j, "md_gto_iso_gpu_extract: out of device memory or buffer overflow"); break; }
            if (nv == 0) {
                isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_OUTPUT, &j->ts);
                j->stage = ISOG_J_DONE;
                break;
            }
            a.out_pos = out->position; a.out_nrm = out->normal; a.out_idx = out->index;
            a.count = nv;
            a.cap = cap_q;      // the kernels take the quad count from counters[1]
            isog_mark(s, "output alloc", nv, cap_q);
            md_gpu_memset(s, a.bucket, 0, sizeof(uint32_t) * 192);
            MD_GPU_LAUNCH(s, isog_k[ISOG_K_BCOUNT], md_gpu_grid_for(isog_k[ISOG_K_BCOUNT], cap_q, 1, 1), a);
            MD_GPU_LAUNCH(s, isog_k[ISOG_K_BSCAN], md_gpu_grid(1, 1, 1), a);
            isog_mark(s, "bucket: count, scan", nv, cap_q);
            MD_GPU_LAUNCH(s, isog_k[ISOG_K_BVERTS], md_gpu_grid_for(isog_k[ISOG_K_BVERTS], nv, 1, 1), a);
            isog_mark(s, "bucket: vertices", nv, 0);
            MD_GPU_LAUNCH(s, isog_k[ISOG_K_BTRIS], md_gpu_grid_for(isog_k[ISOG_K_BTRIS], cap_q, 1, 1), a);
            isog_mark(s, "bucket: triangles", cap_q, 0);
            // the final counts and the buckets; with read_back the mesh too, its triangles up to twice the
            // vertices (a closed surface has about as many quads as vertices; more are read after)
            isog_jfetch_add(j, j->bk, a.bucket, sizeof(j->bk));
            isog_jfetch_add(j, j->ctr, a.counters, sizeof(j->ctr));
            j->idx_copied = 0;
            if (desc->read_back) {
                static int test_short = -1;     // TEST: GTO_ISO_READBACK_SHORT=1 copies no triangles first (the second readback)
                if (test_short < 0) test_short = getenv("GTO_ISO_READBACK_SHORT") != NULL;
                const uint32_t qc = test_short ? 0u : cap_q < 2 * nv + 64 ? cap_q : 2 * nv + 64;
                j->idx_copied = 6 * qc;
                isog_jfetch_add(j, NULL, out->position, sizeof(md_gpu_float4) * nv);
                isog_jfetch_add(j, NULL, out->normal, sizeof(md_gpu_float4) * nv);
                isog_jfetch_add(j, NULL, out->index, sizeof(uint32_t) * j->idx_copied);
            }
            if (!isog_jfetch(j)) { isog_job_fail(j, "md_gto_iso_gpu_extract: out of host-readable memory"); break; }
            if (desc->read_back) {
                j->fetch_keep = true;
                for (int k = 0; k < 3; ++k) j->mesh_off[k] = j->fetch_off[2 + k];
            }
            j->stage = ISOG_J_FINAL_READ;
        } break;

        case ISOG_J_FINAL_READ: {
            const uint32_t* bk = j->bk;
            const uint32_t nv = j->nv;
            isog_jprog(j, 0.95f);
            if (j->kept.gpu) {      // the mesh's host copy: the mesh's from here (freed with it)
                out->host_mem = j->kept.gpu;
                out->host_position = (const float*)((const char*)j->kept.cpu + j->mesh_off[0]);
                out->host_normal = (const float*)((const char*)j->kept.cpu + j->mesh_off[1]);
                out->host_index = (const uint32_t*)((const char*)j->kept.cpu + j->mesh_off[2]);
                j->kept = (md_gpu_mem_t){0};
            }
            isog_mark(s, "readback: counts, buckets", nv, j->ctr[1]);
            if (j->ctr[7]) { isog_job_fail(j, "md_gto_iso_gpu_extract: more than 2048 shells (or 8192 AOs) reach a single point (the D form)"); break; }
            if (j->ctr[3]) { isog_job_fail(j, "md_gto_iso_gpu_extract: out of device memory or buffer overflow"); break; }
            j->nq = j->ctr[1] < j->cap_q ? j->ctr[1] : j->cap_q;
            out->open_edges = j->ctr[2];
            for (uint32_t k = 0; k < F; ++k) {
                for (int sh = 0; sh < 2; ++sh) {
                    const uint32_t i = 2 * k + (uint32_t)sh;
                    out->vertex_offset[k][sh]   = bk[64 + i];
                    out->vertex_count_of[k][sh] = bk[i];
                    out->index_offset[k][sh]    = 6 * bk[96 + i];
                    out->index_count_of[k][sh]  = 6 * bk[32 + i];
                }
            }
            out->vertex_count = nv;
            out->index_count = 6 * j->nq;
            isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_OUTPUT, &j->ts);
            if (!desc->read_back || out->index_count <= j->idx_copied) { j->stage = ISOG_J_DONE; break; }
            // more triangles than the first copy held: the whole mesh again
            const size_t bv = sizeof(float) * 4 * nv, bi = sizeof(uint32_t) * out->index_count;
            isog_jfetch_add(j, NULL, out->position, bv);
            isog_jfetch_add(j, NULL, out->normal, bv);
            isog_jfetch_add(j, NULL, out->index, bi);
            if (!isog_jfetch(j)) { isog_job_fail(j, "md_gto_iso_gpu_extract: out of host-readable memory"); break; }
            j->fetch_keep = true;
            for (int k = 0; k < 3; ++k) j->mesh_off[k] = j->fetch_off[k];
            j->stage = ISOG_J_MESH_READ;
        } break;

        case ISOG_J_MESH_READ:
            md_gpu_free(s, out->host_mem);      // the first copy
            out->host_mem = j->kept.gpu;
            out->host_position = (const float*)((const char*)j->kept.cpu + j->mesh_off[0]);
            out->host_normal = (const float*)((const char*)j->kept.cpu + j->mesh_off[1]);
            out->host_index = (const uint32_t*)((const char*)j->kept.cpu + j->mesh_off[2]);
            j->kept = (md_gpu_mem_t){0};
            j->stage = ISOG_J_DONE;
            break;

        default:
            isog_job_fail(j, "md_gto_iso_gpu_extract: invalid job state");
            break;
        }
    }
#undef a
    if (j->stage == ISOG_J_DONE) {
        isog_jprog(j, 1.0f);
        out->seconds_total = md_tick_to_seconds(md_tick_now() - j->t_start);
    }
}

md_gto_iso_gpu_job_t* md_gto_iso_gpu_job_begin(md_gpu_stream_t s, const md_gto_iso_gpu_set_t* set, const md_gto_iso_gpu_desc_t* desc) {
    const double h = desc ? iso_lattice_float(desc->cell_size) : 0.0;
    if (!s || !set || !desc || !(h > 0.0)) {
        MD_LOG_ERROR("md_gto_iso_gpu_extract: invalid arguments");
        return NULL;
    }
    if (!isog_kernels_ready()) return NULL;
    isog_set_current(s, set);
    const uint32_t F = set->num_fields;
    // root: the union of every field's surface box, as on the CPU (the same lattice for one field)
    double lo[3] = {DBL_MAX, DBL_MAX, DBL_MAX}, hi[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
    isog_field_t fields[MD_GTO_ISO_GPU_MAX_FIELDS];
    memset(fields, 0, sizeof(fields));     // all 16 uploaded: unrolled loops may read past num_fields
    uint32_t root_mask = 0;
    for (uint32_t k = 0; k < F; ++k) {
        if (!(desc->iso[k] > 0.0)) { MD_LOG_ERROR("md_gto_iso_gpu_extract: iso[%u] must be positive", k); return NULL; }
        double bl[3], bh[3];
        iso_surface_box(&set->host[k], desc->iso[k], bl, bh);
        for (int a = 0; a < 3; ++a) { lo[a] = fmin(lo[a], bl[a]); hi[a] = fmax(hi[a], bh[a]); }
        fields[k] = set->fields[k];
        fields[k].iso = (float)desc->iso[k];
        const uint32_t sheets = desc->sheets[k] ? (desc->sheets[k] & 3u) : 3u;
        root_mask |= sheets << (2u * k);
    }
    double origin[3];
    const int L = iso_root_box(lo, hi, 4.0 * h, origin);
    if (L > 8) { MD_LOG_ERROR("md_gto_iso_gpu_extract: lattice too large for cell size %g", h); return NULL; }

    md_gto_iso_gpu_job_t* j = (md_gto_iso_gpu_job_t*)calloc(1, sizeof(md_gto_iso_gpu_job_t));
    if (!j) return NULL;
    j->s = s;
    j->set = set;
    j->desc = *desc;
    j->desc.cell_size = h;
    j->F = F;
    j->L = L;
    memcpy(j->origin, origin, sizeof(origin));
    memcpy(j->fields, fields, sizeof(fields));
    j->out.num_fields = F;
    isog_prof_begin(s);
    j->t_start = j->ts = md_tick_now();
    isog_mark(s, "host: root box", F, (uint32_t)L);
    memcpy(j->out.origin, origin, sizeof(origin));
    j->out.cell_size = h;
    j->out.cells_per_axis = 4u << L;

    isog_args_t* a = &j->a;
    isog_args_set(a, set);
    a->origin_h   = (md_gpu_float4){(float)origin[0], (float)origin[1], (float)origin[2], (float)h};
    a->brick_level = (uint32_t)L;
    a->step_max   = (float)(0.5 * h);
    a->lambda     = 0.5f;
    a->tol        = (float)(pow(10.0, -(double)ISOG_NEWTON_TOL_EXP) * h);
    j->cdform = a->dform != 0;      // the corners in the D form with it
    // the D form's octree bounds the matrix's density, within the factors' slack of what the factored
    // evaluation gives (the slack is in the bound): either evaluation goes with it
    j->odform = set->num_dmat > 0;
    // APPROXIMATE: seed and flood (GTO_ISO_FLOOD): every field (2), or sets of densities whose maxima are at
    // the atoms (md_gto_iso_field_set_maxima_at_atoms) and whose every weight is >= 0 (1): each piece of their
    // surfaces holds an atom, whose brick is a seed
    bool flood = isog_flood == 2;
    if (isog_flood == 1) {
        flood = true;
        for (uint32_t k = 0; k < F && flood; ++k) {
            const md_gto_iso_field_t* fk = set->src[k];
            const iso_density_t* dk = fk ? fk->f.dens : NULL;
            if (!dk || !fk->maxima_at_atoms) { flood = false; break; }
            for (int t = 0; t < dk->num_terms && flood; ++t) flood = dk->terms[t].a == dk->terms[t].b && dk->terms[t].w >= 0.0;
        }
    }
    j->flood = flood;
    a->fields = isog_jalloc(j, sizeof(fields));
    a->counters = isog_jalloc(j, sizeof(uint32_t) * 12);
    a->grid_out = isog_jalloc(j, sizeof(uint32_t) * 4);
    if (!a->fields || !a->counters || !a->grid_out) {
        md_gto_iso_gpu_job_end(j, NULL);
        MD_LOG_ERROR("md_gto_iso_gpu_extract: out of device memory");
        return NULL;
    }
    md_gpu_upload(s, a->fields, fields, sizeof(fields));
    md_gpu_memset(s, a->counters, 0, sizeof(uint32_t) * 12);
    // count_in valid wherever a kernel may read it (llvmpipe loads both sides of a select)
    a->count_in = a->counters;
    j->root_mask = root_mask;
    j->cap = 1u << 18;
    j->stage = ISOG_J_OCTREE;
    isog_job_step(j);        // the GPU starts at once
    return j;
}

md_gto_iso_gpu_job_status_t md_gto_iso_gpu_job_poll(md_gto_iso_gpu_job_t* j, double budget_ms) {
    if (!j) return MD_GTO_ISO_GPU_JOB_FAILED;
    const md_tick_t t0 = md_tick_now();
    for (;;) {
        isog_jmarks_take(j);
        if (j->stage == ISOG_J_DONE) return MD_GTO_ISO_GPU_JOB_DONE;
        if (j->stage == ISOG_J_FAILED) return MD_GTO_ISO_GPU_JOB_FAILED;
        if (j->fetching) {
            if (!md_gpu_sync_is_complete(j->fetch_sync)) {
                if (budget_ms < 0.0) {
                    md_gpu_sync_wait(j->fetch_sync);
                } else if (md_tick_to_seconds(md_tick_now() - t0) * 1e3 >= budget_ms) {
                    return MD_GTO_ISO_GPU_JOB_RUNNING;
                } else {
                    continue;
                }
            }
            isog_jfetch_take(j);
        }
        isog_job_step(j);
    }
}

float md_gto_iso_gpu_job_progress(const md_gto_iso_gpu_job_t* j) {
    return j ? j->progress : 0.0f;
}

bool md_gto_iso_gpu_job_end(md_gto_iso_gpu_job_t* j, md_gto_iso_gpu_mesh_t* out) {
    if (!j) {
        if (out) memset(out, 0, sizeof(*out));
        return false;
    }
    const bool done = j->stage == ISOG_J_DONE;
    if (done && out) {
        *out = j->out;
    } else {
        md_gto_iso_gpu_mesh_free(j->s, &j->out);
        if (out) memset(out, 0, sizeof(*out));
    }
    if (j->fetching && j->fetch_mem.gpu) md_gpu_free(j->s, j->fetch_mem.gpu);   // after the copy in flight
    if (j->kept.gpu) md_gpu_free(j->s, j->kept.gpu);
    for (uint32_t i = 0; i < j->num_alloc; ++i) md_gpu_free(j->s, j->alloc[i]);
    free(j->bl);
    free(j->slot);
    free(j);
    return done && out;
}

bool md_gto_iso_gpu_extract(md_gpu_stream_t s, const md_gto_iso_gpu_set_t* set, const md_gto_iso_gpu_desc_t* desc, md_gto_iso_gpu_mesh_t* out) {
    if (!out) {
        MD_LOG_ERROR("md_gto_iso_gpu_extract: invalid arguments");
        return false;
    }
    memset(out, 0, sizeof(*out));
    md_gto_iso_gpu_job_t* j = md_gto_iso_gpu_job_begin(s, set, desc);
    if (!j) return false;
    j->blocking = true;
    md_gto_iso_gpu_job_poll(j, -1.0);
    return md_gto_iso_gpu_job_end(j, out);
}

void md_gto_iso_gpu_mesh_free(md_gpu_stream_t stream, md_gto_iso_gpu_mesh_t* mesh) {
    if (!mesh) return;
    md_gpu_free(stream, mesh->position);
    md_gpu_free(stream, mesh->normal);
    md_gpu_free(stream, mesh->index);
    if (mesh->host_mem) md_gpu_free(stream, mesh->host_mem);
    memset(mesh, 0, sizeof(*mesh));
}

bool md_gto_iso_gpu_mesh_read(md_gpu_stream_t stream, const md_gto_iso_gpu_mesh_t* mesh, uint32_t field, md_gto_iso_mesh_t* out, struct md_allocator_i* alloc) {
    if (!stream || !mesh || !out || !alloc || field >= mesh->num_fields) return false;
    memset(out, 0, sizeof(*out));
    out->alloc = alloc;
    const uint32_t v0 = mesh->vertex_offset[field][0];
    const uint32_t nv = mesh->vertex_count_of[field][0] + mesh->vertex_count_of[field][1];
    const uint32_t i0 = mesh->index_offset[field][0];
    const uint32_t ni = mesh->index_count_of[field][0] + mesh->index_count_of[field][1];
    out->vertex_count = nv;
    out->index_count = ni;
    out->position = (float*)md_alloc(alloc, sizeof(float) * 3 * (nv + 1));
    out->normal   = (float*)md_alloc(alloc, sizeof(float) * 3 * (nv + 1));
    out->index    = (uint32_t*)md_alloc(alloc, sizeof(uint32_t) * (ni + 1));
    for (int s = 0; s < 2; ++s) {
        out->sheet_vertex_offset[s] = mesh->vertex_offset[field][s] - v0;
        out->sheet_vertex_count[s]  = mesh->vertex_count_of[field][s];
        out->sheet_index_offset[s]  = mesh->index_offset[field][s] - i0;
        out->sheet_index_count[s]   = mesh->index_count_of[field][s];
    }
    if (nv == 0) return true;
    float* p4 = (float*)malloc(sizeof(float) * 4 * nv);
    float* n4 = (float*)malloc(sizeof(float) * 4 * nv);
    bool ok = isog_read(stream, p4, mesh->position + sizeof(float) * 4 * v0, sizeof(float) * 4 * nv)
           && isog_read(stream, n4, mesh->normal + sizeof(float) * 4 * v0, sizeof(float) * 4 * nv)
           && isog_read(stream, out->index, mesh->index + sizeof(uint32_t) * i0, sizeof(uint32_t) * ni);
    for (uint32_t i = 0; i < nv; ++i) {
        for (int k = 0; k < 3; ++k) { out->position[3 * i + k] = p4[4 * i + k]; out->normal[3 * i + k] = n4[4 * i + k]; }
    }
    for (uint32_t i = 0; i < ni; ++i) out->index[i] -= v0;
    free(p4); free(n4);
    return ok;
}

bool md_gto_iso_gpu_mesh_read_all(md_gpu_stream_t stream, const md_gto_iso_gpu_mesh_t* mesh, float* position, float* normal, uint32_t* index) {
    if (!stream || !mesh) return false;
    const size_t bv = sizeof(float) * 4 * mesh->vertex_count;
    const size_t bi = sizeof(uint32_t) * mesh->index_count;
    const size_t bytes = (position ? bv : 0) + (normal ? bv : 0) + (index ? bi : 0);
    if (!bytes) return true;
    md_gpu_mem_t m = md_gpu_malloc(stream, MD_GPU_MEM_HOST_READ, bytes);
    if (!m.gpu) return false;
    size_t off = 0;
    if (position && bv) { md_gpu_copy(stream, m.gpu + off, mesh->position, bv); off += bv; }
    if (normal && bv)   { md_gpu_copy(stream, m.gpu + off, mesh->normal, bv);   off += bv; }
    if (index && bi)    { md_gpu_copy(stream, m.gpu + off, mesh->index, bi);    off += bi; }
    md_gpu_stream_sync(stream);
    off = 0;
    if (position) { memcpy(position, (const char*)m.cpu + off, bv); off += bv; }
    if (normal)   { memcpy(normal,   (const char*)m.cpu + off, bv); off += bv; }
    if (index)    { memcpy(index,    (const char*)m.cpu + off, bi); off += bi; }
    md_gpu_free(stream, m.gpu);
    return true;
}

// ---------------------------------------------------------------------------------------------
// The fields at arbitrary points: the vertex evaluators on the caller's points
// ---------------------------------------------------------------------------------------------

#define ISOG_EVAL_BATCH (1u << 20)      // points per pass: the scratch is ~52 bytes a point

bool md_gto_iso_gpu_evaluate_launch(md_gpu_stream_t s, const md_gto_iso_gpu_set_t* set, const md_gto_iso_gpu_eval_desc_t* desc) {
    if (!s || !set || !desc || !desc->xyz || desc->num_fields == 0 || desc->num_fields > MD_GTO_ISO_GPU_MAX_FIELDS ||
        (!desc->out_value && !desc->out_gradient) || (desc->xyz_stride & 3u)) {
        MD_LOG_ERROR("md_gto_iso_gpu_evaluate: invalid arguments");
        return false;
    }
    for (uint32_t k = 0; k < desc->num_fields; ++k) {
        if (desc->fields[k] >= set->num_fields) {
            MD_LOG_ERROR("md_gto_iso_gpu_evaluate: field %u is not in the set (%u fields)", desc->fields[k], set->num_fields);
            return false;
        }
    }
    if (desc->num_xyz == 0) {
        if (desc->out_failed) md_gpu_memset(s, desc->out_failed, 0, sizeof(uint32_t));
        return true;
    }
    if (!isog_kernels_ready()) return false;
    isog_set_current(s, set);

    isog_args_t a;
    memset(&a, 0, sizeof(a));
    isog_args_set(&a, set);
    const isog_vkernel_t vk = isog_vertex_kernel(set, a.dform != 0u, false);
    if (!isog_have(vk.k) || !isog_have(ISOG_K_PLOAD) || !isog_have(ISOG_K_PSTORE)) return false;

    const uint32_t n = desc->num_xyz;
    const uint32_t B = n < ISOG_EVAL_BATCH ? n : ISOG_EVAL_BATCH;
    const uint32_t stride = desc->xyz_stride ? desc->xyz_stride : 3u * (uint32_t)sizeof(float);
    // the vertices' state; vval the caller's gradients where they are wanted
    md_gpu_mem_t m_fields   = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(isog_field_t) * MD_GTO_ISO_GPU_MAX_FIELDS);
    md_gpu_mem_t m_counters = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(uint32_t) * 12);
    md_gpu_mem_t m_vpos     = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(md_gpu_float4) * B);
    md_gpu_mem_t m_vinfo    = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(md_gpu_uint4) * B);
    md_gpu_mem_t m_vdone    = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(uint32_t) * B);
    md_gpu_mem_t m_vval     = desc->out_gradient ? (md_gpu_mem_t){0} : md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(md_gpu_float4) * B);
    const bool ok = m_fields.gpu && m_counters.gpu && m_vpos.gpu && m_vinfo.gpu && m_vdone.gpu && (desc->out_gradient || m_vval.gpu);
    if (ok) {
        // the fields as the extraction uploads them, all 16 (unrolled loops may read past num_fields); their iso
        // values are not read by an evaluation
        isog_field_t fields[MD_GTO_ISO_GPU_MAX_FIELDS];
        memset(fields, 0, sizeof(fields));
        memcpy(fields, set->fields, sizeof(isog_field_t) * set->num_fields);
        md_gpu_upload(s, m_fields.gpu, fields, sizeof(fields));
        md_gpu_memset(s, m_counters.gpu, 0, sizeof(uint32_t) * 12);
        a.fields   = m_fields.gpu;
        a.counters = m_counters.gpu;
        a.count_in = m_counters.gpu;    // valid wherever a kernel may read it (llvmpipe loads both sides of a select)
        a.vpos     = m_vpos.gpu;
        a.vinfo    = m_vinfo.gpu;
        a.vdone    = m_vdone.gpu;
        a.pts_stride = stride / 4u;
        for (uint32_t off = 0; off < n; off += B) {
            const uint32_t nb = n - off < B ? n - off : B;
            a.count = nb;
            a.pts = desc->xyz + (md_gpu_addr_t)off * stride;
            const md_gpu_grid_t g_pts = md_gpu_grid_for(isog_k[ISOG_K_PLOAD], nb, 1, 1);
            const md_gpu_grid_t g_eval = vk.per_group ? isog_groups((nb + vk.per_group - 1) / vk.per_group) : md_gpu_grid_for(isog_k[vk.k], nb, 1, 1);
            for (uint32_t k = 0; k < desc->num_fields; ++k) {
                const size_t first = (size_t)k * n + off;      // this pass's first output of field k
                a.pts_field = desc->fields[k];
                a.mode = k == 0 ? 1u : 0u;                      // the positions once a pass, then only the field
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_PLOAD], g_pts, a);
                a.mode = 0;
                a.vval = desc->out_gradient ? desc->out_gradient + first * sizeof(md_gpu_float4) : m_vval.gpu;
                md_gpu_memset(s, a.vval, 0xFF, sizeof(md_gpu_float4) * nb);    // NaN where no value is made
                MD_GPU_LAUNCH(s, isog_k[vk.k], g_eval, a);
                if (desc->out_value) {
                    a.pval = desc->out_value + first * sizeof(float);
                    MD_GPU_LAUNCH(s, isog_k[ISOG_K_PSTORE], g_pts, a);
                }
            }
        }
        if (desc->out_failed) md_gpu_copy(s, desc->out_failed, m_counters.gpu + 7 * sizeof(uint32_t), sizeof(uint32_t));
    } else {
        MD_LOG_ERROR("md_gto_iso_gpu_evaluate: out of device memory (%u points)", B);
    }
    // stream ordered: after the work above
    if (m_fields.gpu)   md_gpu_free(s, m_fields.gpu);
    if (m_counters.gpu) md_gpu_free(s, m_counters.gpu);
    if (m_vpos.gpu)     md_gpu_free(s, m_vpos.gpu);
    if (m_vinfo.gpu)    md_gpu_free(s, m_vinfo.gpu);
    if (m_vdone.gpu)    md_gpu_free(s, m_vdone.gpu);
    if (m_vval.gpu)     md_gpu_free(s, m_vval.gpu);
    return ok;
}

// 10 bits of x spread to every third bit
static inline uint32_t isog_morton_spread(uint32_t x) {
    x &= 0x3FFu;
    x = (x | (x << 16)) & 0x030000FFu;
    x = (x | (x << 8))  & 0x0300F00Fu;
    x = (x | (x << 4))  & 0x030C30C3u;
    x = (x | (x << 2))  & 0x09249249u;
    return x;
}

static int isog_cmp_u64(const void* pa, const void* pb) {
    const uint64_t x = *(const uint64_t*)pa, y = *(const uint64_t*)pb;
    return x < y ? -1 : x > y;
}

bool md_gto_iso_gpu_evaluate(md_gpu_stream_t s, const md_gto_iso_gpu_set_t* set, const uint32_t* fields, uint32_t num_fields,
                             const float* xyz, size_t xyz_stride, size_t num_xyz, float* out_value, float* out_gradient) {
    if (!s || !set || !fields || num_fields == 0 || num_fields > MD_GTO_ISO_GPU_MAX_FIELDS || (num_xyz && !xyz) ||
        (!out_value && !out_gradient) || num_xyz > UINT32_MAX) {
        MD_LOG_ERROR("md_gto_iso_gpu_evaluate: invalid arguments");
        return false;
    }
    if (num_xyz == 0) return true;
    const size_t n = num_xyz;
    const size_t st = xyz_stride ? xyz_stride : 3 * sizeof(float);
#define ISOG_P(i) ((const float*)((const char*)xyz + (i) * st))

    // a Morton curve through the points' box: code << 32 | index, sorted
    float lo[3] = {FLT_MAX, FLT_MAX, FLT_MAX}, hi[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
    for (size_t i = 0; i < n; ++i) {
        const float* p = ISOG_P(i);
        for (int c = 0; c < 3; ++c) { lo[c] = MIN(lo[c], p[c]); hi[c] = MAX(hi[c], p[c]); }
    }
    float scl[3];
    for (int c = 0; c < 3; ++c) scl[c] = hi[c] > lo[c] ? 1023.0f / (hi[c] - lo[c]) : 0.0f;
    uint64_t* order = (uint64_t*)malloc(sizeof(uint64_t) * n);
    md_gpu_float4* pts = (md_gpu_float4*)malloc(sizeof(md_gpu_float4) * n);
    if (!order || !pts) {
        free(order); free(pts);
        MD_LOG_ERROR("md_gto_iso_gpu_evaluate: out of memory");
        return false;
    }
    for (size_t i = 0; i < n; ++i) {
        const float* p = ISOG_P(i);
        uint32_t q[3];
        for (int c = 0; c < 3; ++c) {
            const float t = (p[c] - lo[c]) * scl[c];
            q[c] = t > 0.0f ? (t < 1023.0f ? (uint32_t)t : 1023u) : 0u;     // NaN to 0
        }
        const uint32_t code = isog_morton_spread(q[0]) | (isog_morton_spread(q[1]) << 1) | (isog_morton_spread(q[2]) << 2);
        order[i] = ((uint64_t)code << 32) | (uint64_t)i;
    }
    qsort(order, n, sizeof(uint64_t), isog_cmp_u64);
    for (size_t r = 0; r < n; ++r) {
        const float* p = ISOG_P((size_t)(uint32_t)order[r]);
        pts[r] = (md_gpu_float4){p[0], p[1], p[2], 1.0f};
    }
#undef ISOG_P

    // device: the points, the outputs and the count of failures; then all of it to one host-readable buffer
    const size_t bv = out_value ? sizeof(float) * num_fields * n : 0;
    const size_t bg = out_gradient ? sizeof(md_gpu_float4) * num_fields * n : 0;
    md_gpu_mem_t d_pts = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(md_gpu_float4) * n);
    md_gpu_mem_t d_out = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, bv + bg + 16);
    md_gpu_mem_t h_out = md_gpu_malloc(s, MD_GPU_MEM_HOST_READ, bv + bg + 16);
    bool ok = d_pts.gpu && d_out.gpu && h_out.gpu;
    if (ok) {
        md_gpu_upload(s, d_pts.gpu, pts, sizeof(md_gpu_float4) * n);
        md_gto_iso_gpu_eval_desc_t d;
        memset(&d, 0, sizeof(d));
        d.xyz = d_pts.gpu;
        d.xyz_stride = sizeof(md_gpu_float4);
        d.num_xyz = (uint32_t)n;
        d.num_fields = num_fields;
        memcpy(d.fields, fields, sizeof(uint32_t) * num_fields);
        d.out_value    = out_value ? d_out.gpu : 0;
        d.out_gradient = out_gradient ? d_out.gpu + bv : 0;
        d.out_failed   = d_out.gpu + bv + bg;
        ok = md_gto_iso_gpu_evaluate_launch(s, set, &d);
        if (ok) {
            md_gpu_copy(s, h_out.gpu, d_out.gpu, bv + bg + sizeof(uint32_t));
            md_gpu_stream_sync(s);
            const float* hv = (const float*)h_out.cpu;
            const md_gpu_float4* hg = (const md_gpu_float4*)((const char*)h_out.cpu + bv);
            const uint32_t failed = *(const uint32_t*)((const char*)h_out.cpu + bv + bg);
            for (uint32_t k = 0; k < num_fields; ++k) {
                for (size_t r = 0; r < n; ++r) {
                    const size_t i = (size_t)(uint32_t)order[r];
                    if (out_value) out_value[(size_t)k * n + i] = hv[(size_t)k * n + r];
                    if (out_gradient) {
                        const md_gpu_float4 g = hg[(size_t)k * n + r];
                        float* o = out_gradient + 3 * ((size_t)k * n + i);
                        o[0] = g.x; o[1] = g.y; o[2] = g.z;
                    }
                }
            }
            if (failed) {
                MD_LOG_ERROR("md_gto_iso_gpu_evaluate: %u values could not be made (too many shells reach their points)", failed);
                ok = false;
            }
        }
    } else {
        MD_LOG_ERROR("md_gto_iso_gpu_evaluate: out of device memory (%zu points)", n);
    }
    if (d_pts.gpu) md_gpu_free(s, d_pts.gpu);
    if (d_out.gpu) md_gpu_free(s, d_out.gpu);
    if (h_out.gpu) md_gpu_free(s, h_out.gpu);
    free(order);
    free(pts);
    return ok;
}

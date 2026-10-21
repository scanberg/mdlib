// md_gto_iso_gpu.inl - the md_gpu path of md_gto_iso (included by md_gto_iso.c under MD_ENABLE_GPU).
//
// Host side of src/shaders/iso/gto_iso_gpu.slang: builds a field set (basis, coefficient matrix C,
// terms per field), then issues the pipeline: octree level by level, bricks, corner values, Surface
// Nets, Newton projection and relaxation, per-field regrouping. Three waits for the GPU: the brick
// list (to size the corner storage), the vertex count (to size the mesh buffers), the final counts.

#include <core/md_intrinsics.h>
#include <core/md_gpu.h>
#include <gto_iso_gpu_shaders.inl>

// BENCHMARK (temporary, branch gto-iso-gpu-opt): which round of changes this is, in every profile report
#define ISOG_REVISION "r25"

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
    uint32_t canon, diag;
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
} isog_args_t;

enum {
    ISOG_K_GRID, ISOG_K_OCTREE, ISOG_K_HASH, ISOG_K_CORNERS, ISOG_K_EVAL, ISOG_K_CELLS, ISOG_K_QUADS,
    ISOG_K_NEWTON, ISOG_K_RELAX, ISOG_K_BCOUNT, ISOG_K_BSCAN, ISOG_K_BVERTS, ISOG_K_BTRIS, ISOG_K_EVAL1, ISOG_K_CORNERS1,
    ISOG_K_CORNERS_F1, ISOG_K_CORNERS_F8, ISOG_K_CORNERS_F16, ISOG_K_CORNERS_F32, ISOG_K_CORNERS_F64,
    ISOG_K_EVAL_F8, ISOG_K_EVAL_F16, ISOG_K_EVAL_F32, ISOG_K_EVAL_F64, ISOG_K_CORNERS_FILL,
    ISOG_K_POINTS, ISOG_K_CLIST_F1, ISOG_K_CLIST_F8, ISOG_K_CLIST_F16, ISOG_K_CLIST_F32, ISOG_K_CLIST_F64,
    ISOG_K_OCTREE_SMALL, ISOG_K_GCORN_F16, ISOG_K_GCORN_F32, ISOG_K_GCORN_F64, ISOG_K_GEVAL_F16, ISOG_K_GEVAL_F32, ISOG_K_GEVAL_F64,
    ISOG_K_DCORN, ISOG_K_DEVAL, ISOG_K_ODFORM, ISOG_K_EVAL_LIST1,
    ISOG_K_GCORN_F64_P1, ISOG_K_GCORN_F64_P2, ISOG_K_GEVAL_F64_P1, ISOG_K_GEVAL_F64_P2, ISOG_K_CLIST_F1V4, ISOG_K_DCORN2, ISOG_K_DEVAL2,
    ISOG_K_DCORN3, ISOG_K_DEVAL3, ISOG_K_SCAN,
    ISOG_K_GCORN_W128, ISOG_K_GCORN_W192, ISOG_K_GCORN_W256, ISOG_K_GEVAL_W128, ISOG_K_GEVAL_W192, ISOG_K_GEVAL_W256,
    ISOG_K_ODFORM8,
    ISOG_K_GCORN_X128, ISOG_K_GCORN_X192, ISOG_K_GCORN_X256, ISOG_K_GEVAL_X128, ISOG_K_GEVAL_X192, ISOG_K_GEVAL_X256,
    ISOG_K_GCORN_Y128, ISOG_K_GCORN_Y192, ISOG_K_GCORN_Y256, ISOG_K_GEVAL_Y128, ISOG_K_GEVAL_Y192, ISOG_K_GEVAL_Y256,
    ISOG_K_FSEED, ISOG_K_FFACE, ISOG_K_FFINISH,
    ISOG_K_COUNT
};

static md_gpu_kernel_t isog_k[ISOG_K_COUNT];
static const char*     isog_kname[ISOG_K_COUNT];   // their labels

static md_gpu_device_t isog_device;
static bool             isog_failed;

// EXPERIMENT (temporary, branch gto-iso-gpu-opt): switches between the old and the new forms, also
// from the environment (read at initialization) so one build measures both:
//   GTO_ISO_RELAX_FIRST=0   project, then 3 x (relax with the field normals, project)
//                       1   relax on the lattice surface (trilinear normals), then project (default)
//   GTO_ISO_CORNERS_V2=0    iso_corners (GEMM tiles) and iso_corners_orbital
//                      1    the register form, iso_corners_f* (default)
//   GTO_ISO_EVAL_V2=0       iso_eval (GEMM tiles) at the vertices of densities
//                   1       the register form, iso_eval_f* (default)
static int isog_relax_first = 1;
static int isog_corners_v2 = 1;
static int isog_eval_v2 = 1;
//   GTO_ISO_NEWTON_PREDICT=1 a Newton step whose next error is predicted far below the tolerance ends
//                            the vertex without an evaluation to check it (default 1)
static int isog_newton_predict = 1;
//   GTO_ISO_SKIP_BRICK_TEST=1 the octree's last level (the bricks) is not bounded: every child of a kept
//                             node is a brick (default 0)
static int isog_skip_brick_test = 0;
//   GTO_ISO_SHARED_CORNERS=1  a brick's high-face corners taken from the neighbour that owns them when it
//                             is present, instead of evaluated twice (default 1; with the register form)
static int isog_shared_corners = 1;
//   GTO_ISO_EVAL_EPS=10|8|6   the AO cutoff 1e-N of the evaluation (and of the octree's culling, whose
//                             bound adds the dropped part): the volume path uses 1e-6 (default 6)
static int isog_eval_eps_exp = 6;
//   GTO_ISO_CORNERS_LIST=1    corners from a list of the points to evaluate, 256 consecutive points (2-4
//                             neighbouring bricks) per group (default 1; with the register form)
static int isog_corners_list = 1;
//   GTO_ISO_CANON_TERMS=1     the register-form epilogues apply the fields' terms from per-column weights
//                             (w a, w a^2, w a b for adjacent pairs) held by every thread, when every
//                             term has one of those forms (default 1); else from the terms list
static int isog_canon_terms = 1;
//   GTO_ISO_GEMM=1            corners (point list) and density vertices as a tiled GEMM: 128 rows per group,
//                             the AO values of a chunk in groupshared memory, 4 rows x F / 4 columns per
//                             thread (default 1)
static int isog_gemm = 1;
//   GTO_ISO_OCTREE_SMALL=1    sets of at most 4 functions: the octree with a thread per shell, its sums in
//                             registers (default 1)
static int isog_octree_small = 1;
//   GTO_ISO_DFORM=1           sets of densities: corners (point list) and vertices as phi^T D phi over the AOs
//                             reaching each group, D from the fields' terms (default 0: on C60 at A ~ r it is
//                             slower than the factored GEMM; sets beyond the factored kernels' 1024 functions
//                             use it regardless)
static int isog_dform = 0;
//   GTO_ISO_DFORM_OCTREE=1    sets of densities: the octree bounds phi^T D phi (centred form over the AOs a node
//                             reaches; no factors, so no limit on the rank), whichever evaluation (default 1)
static int isog_dform_octree = 1;
//   GTO_ISO_COEF_SCREEN=1     sets of orbitals: a shell is dropped where its largest coefficient times its AOs is
//                             below eps (as the volume path screens), not its AOs alone; the bound's tail per
//                             function becomes 2 eps n_AO (default 1)
static int isog_coef_screen = 1;
//   GTO_ISO_EVAL_LIST=1       one orbital: the vertices' values from the shells listed for 128 vertices at a
//                             time (iso_eval_list_f1), not every shell per vertex (iso_eval_orbital) (default 1)
static int isog_eval_list = 1;
//   GTO_ISO_GEMM_T32=1        the GEMM kernels in tiles of 32 functions for sets of more (fewer registers and
//                             less groupshared memory per group; the AO values made twice as often) (default 0)
static int isog_gemm_t32 = 0;
//   GTO_ISO_GEMM_PIPE=0|1|2   the F = 64 GEMM kernels' loop over a chunk's rows: 0 left to the compiler, 1 not
//                             unrolled, 2 unrolled over 16 rows (padded with zeros) (default 0; what the
//                             driver's register allocation makes of each is the question)
static int isog_gemm_pipe = 0;
//   GTO_ISO_CORNERS_V4=1      one function: four corner points per thread (iso_corners_list_f1v4), a shell's
//                             data read once for the four (default 1)
static int isog_corners_v4 = 1;
//   GTO_ISO_NEWTON_TOL=4|3|2  a vertex has converged within 1e-N h of the surface (default 2: the facets' own
//                             deviation from the surface is about 3e-2 h rms, and on amide 1e-2 h moves it by
//                             at most 1.3% of that; on C60 at h 0.2 it saves the density's second evaluation)
static int isog_newton_tol_exp = 2;
//   GTO_ISO_DFORM_KERNEL=0|1|2  the D form's kernels: 0 those of r9 (iso_dcorners, iso_deval); 1 D' upper
//                             triangular by chunks so Phi . T needs no AO values made again, the chunks' shells
//                             from a table and staged a chunk ahead, a wave's zero chunks skipped, tiles of 64
//                             columns (iso_dcorners2, iso_deval2); 2 the same in tiles of 128 columns, groups of
//                             256 threads (iso_dcorners3, iso_deval3): half the chunks per row (default 2: on
//                             the GTX 1060, C60 density, corners 365 ms against 533 (r9) and 820 (1))
static int isog_dform_kernel = 2;
//   GTO_ISO_ORDERED=1         the bricks sorted by Morton code (siblings stay together), the corner point list
//                             and the vertices in that order (counts, a prefix sum, then offsets) instead of the
//                             order the groups' atomics happen in: consecutive points and vertices close
//                             together, so a group's box (and the shells it lists) stays small (default 1)
static int isog_ordered = 1;
//   GTO_ISO_DIAG=0..5         DIAGNOSTIC, the corner GEMM kernels only (iso_gcorners_*, iso_dcorners2/3): 1 the
//                             AO values not made (a constant), 2 no multiply-adds; iso_gcorners_*: 3 the shell
//                             lists only (no chunks), 4 no rows of C staged, 5 no shells staged and no AO values.
//                             The values are wrong (no surface): only the corner kernel's time means anything
//                             (default 0)
static int isog_diag = 0;
//   GTO_ISO_GEMM_WIDE=1       canonical sets of more than 64 functions: the GEMM in passes of 128, 192 or 256
//                             columns, 64 rows per group of 256 threads (iso_gcorners_w*, iso_geval_w*), so the
//                             shells are listed and staged and the AO values made once per pass instead of once
//                             per tile of 64 (default 1; TEST: 2 for sets of more than 16 functions)
static int isog_gemm_wide = 1;
//   GTO_ISO_OCTREE8=1         the D form's octree below the root with a node's 8 children in one group
//                             (iso_octree_dform8): the shells reaching the parent listed once, each element of
//                             D loaded once for the 8 instead of once per child; the same sums. Off by default
//                             since r22: on the M4 Pro (mm density, r21) the levels of few nodes ran a group
//                             where there were 8 (octree 15.1 ms against 11.4), the bricks' level 6% faster
static int isog_octree8 = 0;
//   GTO_ISO_OCTREE_START=2    the first level the octree tests: the levels above it keep every node (the root
//                             box holds the surfaces, so levels 0 and 1 never dropped one in the series). Their
//                             few nodes reach every AO and take a group each, one after the other: on the M4
//                             Pro the D form's levels 0 and 1 took 2.9 ms (default 2)
static int isog_octree_start = 2;
//   GTO_ISO_DFORM_CORNERS=1   sets of densities: the corners in the D form (the D-form kernel chosen) and the
//                             vertices as without it (the canonical GEMM where there are canonical terms)
//                             (default 0)
static int isog_dform_corners = 0;
//   GTO_ISO_FLOOD=0|1|2       APPROXIMATE: the corners from seeds and a flood instead of every brick the octree
//                             keeps. The bricks' corners are evaluated first; one whose corners straddle an iso
//                             value is a seed, evaluated in full; then, round by round, a brick across a face of
//                             an evaluated one with a sign change on that face. Only bricks with surface cells
//                             are evaluated, but a surface piece no seed reaches (smaller than about 4 cells, or
//                             between the brick corners) is missed. 1: densities with every weight >= 0 (total,
//                             alpha, beta); 2: every field (to see what it misses); 0 off (default 0)
static int isog_flood = 0;
//   GTO_ISO_WIDE_KERNEL=0|1|2 the wide GEMM's kernels: 0 those of r16 (iso_gcorners_w*, iso_geval_w*); 1 the
//                             group's shells listed once (all, in rounds of 256) and staged as they are listed,
//                             chunks of 16 consecutive AOs of the list, so a chunk stages only its rows of C
//                             (iso_gcorners_x*, iso_geval_x*); 2 the same with the next chunk's rows of C loaded
//                             into registers during this chunk's AO values and products (iso_gcorners_y*,
//                             iso_geval_y*: F / 16 more registers) (default 1)
static int isog_wide_kernel = 1;

// The wide kernel k (one of the six ISOG_K_GCORN_W* / ISOG_K_GEVAL_W*) in the form chosen, where the
// device could make it
static int isog_wide_form(int k) {
    const int base = isog_wide_kernel == 2 ? ISOG_K_GCORN_Y128 : isog_wide_kernel == 1 ? ISOG_K_GCORN_X128 : -1;
    if (base < 0) return k;
    const int alt = base + (k - ISOG_K_GCORN_W128);
    return isog_k[alt] ? alt : k;
}

// The kernels are made on first use, so initializing costs nothing until a mesh is extracted
void md_gto_iso_gpu_initialize(md_gpu_device_t device) {
    isog_device = device;
    isog_failed = false;
    const char* e = getenv("GTO_ISO_RELAX_FIRST");
    if (e && e[0]) isog_relax_first = atoi(e);
    e = getenv("GTO_ISO_CORNERS_V2");
    if (e && e[0]) isog_corners_v2 = atoi(e);
    e = getenv("GTO_ISO_EVAL_V2");
    if (e && e[0]) isog_eval_v2 = atoi(e);
    e = getenv("GTO_ISO_NEWTON_PREDICT");
    if (e && e[0]) isog_newton_predict = atoi(e);
    e = getenv("GTO_ISO_SKIP_BRICK_TEST");
    if (e && e[0]) isog_skip_brick_test = atoi(e);
    e = getenv("GTO_ISO_SHARED_CORNERS");
    if (e && e[0]) isog_shared_corners = atoi(e);
    e = getenv("GTO_ISO_EVAL_EPS");
    if (e && e[0]) isog_eval_eps_exp = atoi(e);
    e = getenv("GTO_ISO_CORNERS_LIST");
    if (e && e[0]) isog_corners_list = atoi(e);
    e = getenv("GTO_ISO_CANON_TERMS");
    if (e && e[0]) isog_canon_terms = atoi(e);
    e = getenv("GTO_ISO_GEMM");
    if (e && e[0]) isog_gemm = atoi(e);
    e = getenv("GTO_ISO_OCTREE_SMALL");
    if (e && e[0]) isog_octree_small = atoi(e);
    e = getenv("GTO_ISO_DFORM");
    if (e && e[0]) isog_dform = atoi(e);
    e = getenv("GTO_ISO_DFORM_OCTREE");
    if (e && e[0]) isog_dform_octree = atoi(e);
    e = getenv("GTO_ISO_COEF_SCREEN");
    if (e && e[0]) isog_coef_screen = atoi(e);
    e = getenv("GTO_ISO_EVAL_LIST");
    if (e && e[0]) isog_eval_list = atoi(e);
    e = getenv("GTO_ISO_GEMM_T32");
    if (e && e[0]) isog_gemm_t32 = atoi(e);
    e = getenv("GTO_ISO_GEMM_PIPE");
    if (e && e[0]) isog_gemm_pipe = atoi(e);
    e = getenv("GTO_ISO_CORNERS_V4");
    if (e && e[0]) isog_corners_v4 = atoi(e);
    e = getenv("GTO_ISO_NEWTON_TOL");
    if (e && e[0]) isog_newton_tol_exp = atoi(e);
    if (isog_newton_tol_exp < 1 || isog_newton_tol_exp > 6) isog_newton_tol_exp = 2;
    e = getenv("GTO_ISO_ORDERED");
    if (e && e[0]) isog_ordered = atoi(e);
    e = getenv("GTO_ISO_DIAG");
    if (e && e[0]) isog_diag = atoi(e);
    e = getenv("GTO_ISO_GEMM_WIDE");
    if (e && e[0]) isog_gemm_wide = atoi(e);
    e = getenv("GTO_ISO_OCTREE8");
    if (e && e[0]) isog_octree8 = atoi(e);
    e = getenv("GTO_ISO_OCTREE_START");
    if (e && e[0]) isog_octree_start = atoi(e);
    if (isog_octree_start < 0 || isog_octree_start > 3) isog_octree_start = 2;
    e = getenv("GTO_ISO_WIDE_KERNEL");
    if (e && e[0]) isog_wide_kernel = atoi(e);
    if (isog_wide_kernel < 0 || isog_wide_kernel > 2) isog_wide_kernel = 1;
    e = getenv("GTO_ISO_DFORM_KERNEL");
    if (e && e[0]) isog_dform_kernel = atoi(e);
    if (isog_dform_kernel < 0 || isog_dform_kernel > 2) isog_dform_kernel = 2;
    e = getenv("GTO_ISO_DFORM_CORNERS");
    if (e && e[0]) isog_dform_corners = atoi(e);
    e = getenv("GTO_ISO_FLOOD");
    if (e && e[0]) isog_flood = atoi(e);
    if (isog_flood < 0 || isog_flood > 2) isog_flood = 0;
}

// Every kernel is made at once. One a device cannot make (more groupshared memory than it allows, say)
// fails the extraction only where it is chosen (the D form falls back to its older kernels, the wide
// GEMM to tiles of 64); the kernels every extraction needs fail it always.
static bool isog_ready;

static bool isog_kernel_core(int i) {
    return i == ISOG_K_GRID || i == ISOG_K_HASH || i == ISOG_K_CELLS || i == ISOG_K_QUADS || i == ISOG_K_NEWTON ||
           i == ISOG_K_RELAX || i == ISOG_K_BCOUNT || i == ISOG_K_BSCAN || i == ISOG_K_BVERTS || i == ISOG_K_BTRIS ||
           i == ISOG_K_POINTS || i == ISOG_K_CORNERS_FILL || i == ISOG_K_SCAN;
}

// Kernel k, or false (with a message) when this device could not make it
static bool isog_have(int k) {
    if (k >= 0 && k < ISOG_K_COUNT && isog_k[k]) return true;
    MD_LOG_ERROR("md_gto_iso_gpu_extract: kernel %s is not available on this device", (k >= 0 && k < ISOG_K_COUNT && isog_kname[k]) ? isog_kname[k] : "?");
    return false;
}

static bool isog_kernels_ready(void) {
    if (isog_ready) return true;
    if (isog_failed) return false;
    if (!isog_device) { MD_LOG_ERROR("md_gto_iso_gpu_extract: md_gto_iso_gpu_initialize not called"); return false; }
    const md_gpu_kernel_desc_t d[ISOG_K_COUNT] = {
        md_shader_gto_iso_gpu_iso_grid_kernel(),
        md_shader_gto_iso_gpu_iso_octree_kernel(),
        md_shader_gto_iso_gpu_iso_hash_kernel(),
        md_shader_gto_iso_gpu_iso_corners_kernel(),
        md_shader_gto_iso_gpu_iso_eval_kernel(),
        md_shader_gto_iso_gpu_iso_cells_kernel(),
        md_shader_gto_iso_gpu_iso_quads_kernel(),
        md_shader_gto_iso_gpu_iso_newton_kernel(),
        md_shader_gto_iso_gpu_iso_relax_kernel(),
        md_shader_gto_iso_gpu_iso_bucket_count_kernel(),
        md_shader_gto_iso_gpu_iso_bucket_scan_kernel(),
        md_shader_gto_iso_gpu_iso_bucket_vertices_kernel(),
        md_shader_gto_iso_gpu_iso_bucket_triangles_kernel(),
        md_shader_gto_iso_gpu_iso_eval_orbital_kernel(),
        md_shader_gto_iso_gpu_iso_corners_orbital_kernel(),
        md_shader_gto_iso_gpu_iso_corners_f1_kernel(),
        md_shader_gto_iso_gpu_iso_corners_f8_kernel(),
        md_shader_gto_iso_gpu_iso_corners_f16_kernel(),
        md_shader_gto_iso_gpu_iso_corners_f32_kernel(),
        md_shader_gto_iso_gpu_iso_corners_f64_kernel(),
        md_shader_gto_iso_gpu_iso_eval_f8_kernel(),
        md_shader_gto_iso_gpu_iso_eval_f16_kernel(),
        md_shader_gto_iso_gpu_iso_eval_f32_kernel(),
        md_shader_gto_iso_gpu_iso_eval_f64_kernel(),
        md_shader_gto_iso_gpu_iso_corners_fill_kernel(),
        md_shader_gto_iso_gpu_iso_points_kernel(),
        md_shader_gto_iso_gpu_iso_corners_list_f1_kernel(),
        md_shader_gto_iso_gpu_iso_corners_list_f8_kernel(),
        md_shader_gto_iso_gpu_iso_corners_list_f16_kernel(),
        md_shader_gto_iso_gpu_iso_corners_list_f32_kernel(),
        md_shader_gto_iso_gpu_iso_corners_list_f64_kernel(),
        md_shader_gto_iso_gpu_iso_octree_small_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_f16_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_f32_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_f64_kernel(),
        md_shader_gto_iso_gpu_iso_geval_f16_kernel(),
        md_shader_gto_iso_gpu_iso_geval_f32_kernel(),
        md_shader_gto_iso_gpu_iso_geval_f64_kernel(),
        md_shader_gto_iso_gpu_iso_dcorners_kernel(),
        md_shader_gto_iso_gpu_iso_deval_kernel(),
        md_shader_gto_iso_gpu_iso_octree_dform_kernel(),
        md_shader_gto_iso_gpu_iso_eval_list_f1_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_f64_p1_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_f64_p2_kernel(),
        md_shader_gto_iso_gpu_iso_geval_f64_p1_kernel(),
        md_shader_gto_iso_gpu_iso_geval_f64_p2_kernel(),
        md_shader_gto_iso_gpu_iso_corners_list_f1v4_kernel(),
        md_shader_gto_iso_gpu_iso_dcorners2_kernel(),
        md_shader_gto_iso_gpu_iso_deval2_kernel(),
        md_shader_gto_iso_gpu_iso_dcorners3_kernel(),
        md_shader_gto_iso_gpu_iso_deval3_kernel(),
        md_shader_gto_iso_gpu_iso_scan_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_w128_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_w192_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_w256_kernel(),
        md_shader_gto_iso_gpu_iso_geval_w128_kernel(),
        md_shader_gto_iso_gpu_iso_geval_w192_kernel(),
        md_shader_gto_iso_gpu_iso_geval_w256_kernel(),
        md_shader_gto_iso_gpu_iso_octree_dform8_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_x128_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_x192_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_x256_kernel(),
        md_shader_gto_iso_gpu_iso_geval_x128_kernel(),
        md_shader_gto_iso_gpu_iso_geval_x192_kernel(),
        md_shader_gto_iso_gpu_iso_geval_x256_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_y128_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_y192_kernel(),
        md_shader_gto_iso_gpu_iso_gcorners_y256_kernel(),
        md_shader_gto_iso_gpu_iso_geval_y128_kernel(),
        md_shader_gto_iso_gpu_iso_geval_y192_kernel(),
        md_shader_gto_iso_gpu_iso_geval_y256_kernel(),
        md_shader_gto_iso_gpu_iso_flood_seed_kernel(),
        md_shader_gto_iso_gpu_iso_flood_face_kernel(),
        md_shader_gto_iso_gpu_iso_flood_finish_kernel(),
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
}

// ---------------------------------------------------------------------------------------------
// The field set
// ---------------------------------------------------------------------------------------------

#define ISOG_FN_CAP 1024     // O_FN_CAP in the shader

struct md_gto_iso_gpu_set_t {
    md_gpu_mem_t shells, prims, crit, ao_nrm, ao_ijk, C, terms, fn_tail, wts, dmat, fdmat, drs, dprop;
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
    int           screen;           // the radii are coefficient-screened (orbitals with isog_coef_screen)
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

md_gto_iso_gpu_set_t* md_gto_iso_gpu_set_create(md_gpu_stream_t stream, const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
                                                const md_gto_iso_field_t* const* fields, uint32_t num_fields, struct md_allocator_i* alloc) {
    if (!stream || !basis || !atom_xyz || !fields || !alloc || num_fields == 0 || num_fields > MD_GTO_ISO_GPU_MAX_FIELDS) {
        MD_LOG_ERROR("md_gto_iso_gpu_set_create: invalid arguments (1 to %d fields)", MD_GTO_ISO_GPU_MAX_FIELDS);
        return NULL;
    }
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
    if (dform) {
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
            } else {
                // the terms' symmetrized products
                memset(Dd, 0, sizeof(double) * (size_t)nao * nao);
                for (int t = 0; t < d->num_terms; ++t) {
                    const iso_term_t* tm = &d->terms[t];
                    const double w = tm->a == tm->b ? tm->w : 0.5 * tm->w;
                    for (uint32_t i = 0; i < nao; ++i) {
                        const double ca = d->C[(size_t)i * d->num_fn + tm->a], cb = d->C[(size_t)i * d->num_fn + tm->b];
                        if (ca == 0.0 && cb == 0.0) continue;
                        double* row = Dd + (size_t)i * nao;
                        if (tm->a == tm->b) {
                            const double wa = w * ca;
                            for (uint32_t j = 0; j < nao; ++j) row[j] += wa * d->C[(size_t)j * d->num_fn + tm->a];
                        } else {
                            for (uint32_t j = 0; j < nao; ++j) row[j] += w * (ca * d->C[(size_t)j * d->num_fn + tm->b] + cb * d->C[(size_t)j * d->num_fn + tm->a]);
                        }
                    }
                }
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
        if (!Dm || !Rs) dform = false;
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
        if (num_dmat) {
            ISOG_UP(dmat, Dm, sizeof(float) * (size_t)num_dmat * nao * nao);
            ISOG_UP(fdmat, fdm, sizeof(fdm));
            ISOG_UP(drs, Rs, sizeof(float) * (size_t)num_dmat * nao);
            set->dprop = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(md_gpu_float4) * MD_GTO_ISO_GPU_MAX_FIELDS);
            ok &= set->dprop.gpu != 0;
        }
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

// Every third bit of v (21 bits): one axis of a Morton code
static uint64_t isog_morton_spread(uint32_t v) {
    uint64_t x = v & 0x1FFFFFu;
    x = (x | x << 32) & 0x1F00000000FFFFull;
    x = (x | x << 16) & 0x1F0000FF0000FFull;
    x = (x | x << 8)  & 0x100F00F00F00F00Full;
    x = (x | x << 4)  & 0x10C30C30C30C30C3ull;
    x = (x | x << 2)  & 0x1249249249249249ull;
    return x;
}

// The bricks in Morton order (a stable radix sort): a parent's eight children stay together, and
// bricks close in the list are close in space
static bool isog_sort_bricks(md_gpu_uint4* bl, uint32_t nb) {
    if (nb < 2) return true;
    uint64_t* key = (uint64_t*)malloc(sizeof(uint64_t) * 2 * nb);
    uint32_t* idx = (uint32_t*)malloc(sizeof(uint32_t) * 2 * nb);
    md_gpu_uint4* tmp = (md_gpu_uint4*)malloc(sizeof(md_gpu_uint4) * nb);
    bool ok = key && idx && tmp;
    if (ok) {
        uint64_t kmax = 0;
        for (uint32_t i = 0; i < nb; ++i) {
            key[i] = isog_morton_spread(bl[i].x) | (isog_morton_spread(bl[i].y) << 1) | (isog_morton_spread(bl[i].z) << 2);
            idx[i] = i;
            kmax |= key[i];
        }
        uint64_t* k0 = key; uint64_t* k1 = key + nb;
        uint32_t* i0 = idx; uint32_t* i1 = idx + nb;
        for (int sh = 0; sh < 64 && (kmax >> sh) != 0; sh += 8) {
            uint32_t cnt[257] = {0};
            for (uint32_t i = 0; i < nb; ++i) cnt[((k0[i] >> sh) & 255u) + 1]++;
            for (int d = 0; d < 256; ++d) cnt[d + 1] += cnt[d];
            for (uint32_t i = 0; i < nb; ++i) {
                const uint32_t at = cnt[(k0[i] >> sh) & 255u]++;
                k1[at] = k0[i];
                i1[at] = i0[i];
            }
            uint64_t* kt = k0; k0 = k1; k1 = kt;
            uint32_t* it = i0; i0 = i1; i1 = it;
        }
        for (uint32_t i = 0; i < nb; ++i) tmp[i] = bl[i0[i]];
        memcpy(bl, tmp, sizeof(md_gpu_uint4) * nb);
    }
    free(key);
    free(idx);
    free(tmp);
    return ok;
}

static md_gpu_addr_t isog_temp(md_gpu_stream_t s, size_t bytes) {
    return md_gpu_temp_alloc(s, MD_GPU_MEM_DEVICE, bytes ? bytes : 16).gpu;
}

static md_gpu_grid_t isog_groups(uint32_t n) {     // one group per item, split over y
    if (n == 0) return md_gpu_grid(0, 0, 0);
    return md_gpu_grid(n < 65535u ? n : 65535u, (n + 65534u) / 65535u, 1);
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
    out->relax_first = isog_relax_first != 0;
    out->corners_v2 = isog_corners_v2 != 0;
    out->eval_v2 = isog_eval_v2 != 0;
    out->newton_predict = isog_newton_predict != 0;
    out->skip_brick_test = isog_skip_brick_test != 0;
    out->shared_corners = isog_shared_corners != 0;
    out->eval_eps_exp = isog_eval_eps_exp;
    out->corners_list = isog_corners_list != 0;
    out->canon_terms = isog_canon_terms != 0;
    out->gemm = isog_gemm != 0;
    out->octree_small = isog_octree_small != 0;
    out->dform = isog_dform != 0;
    out->dform_octree = isog_dform_octree != 0;
    out->coef_screen = isog_coef_screen != 0;
    out->eval_list = isog_eval_list != 0;
    out->gemm_t32 = isog_gemm_t32 != 0;
    out->gemm_pipe = isog_gemm_pipe;
    out->corners_v4 = isog_corners_v4 != 0;
    out->newton_tol_exp = isog_newton_tol_exp;
    out->dform_kernel = isog_dform_kernel;
    out->ordered = isog_ordered != 0;
    out->diag = isog_diag;
    out->gemm_wide = isog_gemm_wide != 0;
    out->octree8 = isog_octree8 != 0;
    out->wide_kernel = isog_wide_kernel;
    out->octree_start = isog_octree_start;
    out->dform_corners = isog_dform_corners != 0;
    out->flood = isog_flood;
}

void md_gto_iso_gpu_switches_set(const md_gto_iso_gpu_switches_t* in) {
    if (!in) return;
    isog_relax_first = in->relax_first;
    isog_corners_v2 = in->corners_v2;
    isog_eval_v2 = in->eval_v2;
    isog_newton_predict = in->newton_predict;
    isog_skip_brick_test = in->skip_brick_test;
    isog_shared_corners = in->shared_corners;
    if (in->eval_eps_exp >= 4 && in->eval_eps_exp <= 12) isog_eval_eps_exp = in->eval_eps_exp;
    isog_corners_list = in->corners_list;
    isog_canon_terms = in->canon_terms;
    isog_gemm = in->gemm;
    isog_octree_small = in->octree_small;
    isog_dform = in->dform;
    isog_dform_octree = in->dform_octree;
    isog_coef_screen = in->coef_screen;
    isog_eval_list = in->eval_list;
    isog_gemm_t32 = in->gemm_t32;
    if (in->gemm_pipe >= 0 && in->gemm_pipe <= 2) isog_gemm_pipe = in->gemm_pipe;
    isog_corners_v4 = in->corners_v4;
    if (in->newton_tol_exp >= 1 && in->newton_tol_exp <= 6) isog_newton_tol_exp = in->newton_tol_exp;
    if (in->dform_kernel >= 0 && in->dform_kernel <= 2) isog_dform_kernel = in->dform_kernel;
    isog_dform_corners = in->dform_corners;
    if (in->flood >= 0 && in->flood <= 2) isog_flood = in->flood;
    isog_ordered = in->ordered;
    isog_diag = (in->diag >= 0 && in->diag <= 5) ? in->diag : 0;
    isog_gemm_wide = in->gemm_wide ? (isog_gemm_wide ? isog_gemm_wide : 1) : 0;
    isog_octree8 = in->octree8;
    if (in->wide_kernel >= 0 && in->wide_kernel <= 2) isog_wide_kernel = in->wide_kernel;
    if (in->octree_start >= 0 && in->octree_start <= 3) isog_octree_start = in->octree_start;
}

const char* md_gto_iso_gpu_profile_variant(void) {
    static char buf[512];
    snprintf(buf, sizeof(buf), "%s | relax_first %d corners_v2 %d corners_list %d shared_corners %d eval_v2 %d canon_terms %d gemm %d octree_small %d dform %d dform_kernel %d ordered %d gemm_wide %d wide_kernel %d octree8 %d octree_start %d dform_corners %d flood %d dform_octree %d coef_screen %d eval_list %d gemm_t32 %d gemm_pipe %d corners_v4 %d newton_predict %d newton_tol 1e-%d h skip_brick_test %d eval_eps 1e-%d%s",
             ISOG_REVISION, isog_relax_first, isog_corners_v2, isog_corners_list, isog_shared_corners, isog_eval_v2, isog_canon_terms,
             isog_gemm, isog_octree_small, isog_dform, isog_dform_kernel, isog_ordered, isog_gemm_wide, isog_wide_kernel, isog_octree8, isog_octree_start, isog_dform_corners, isog_flood, isog_dform_octree, isog_coef_screen, isog_eval_list, isog_gemm_t32, isog_gemm_pipe, isog_corners_v4, isog_newton_predict, isog_newton_tol_exp, isog_skip_brick_test, isog_eval_eps_exp,
             isog_diag == 1 ? " DIAG: no AO values" : isog_diag == 2 ? " DIAG: no multiply-adds" : isog_diag == 3 ? " DIAG: lists only"
             : isog_diag == 4 ? " DIAG: no C rows staged" : isog_diag == 5 ? " DIAG: no shells staged, no AO values" : "");
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
        if (strcmp(name[k], "octree level") == 0 && items[k] > 0 && (p->kernels_used & ((1ull << ISOG_K_ODFORM) | (1ull << ISOG_K_ODFORM8)))) {
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

bool md_gto_iso_gpu_extract(md_gpu_stream_t s, const md_gto_iso_gpu_set_t* set, const md_gto_iso_gpu_desc_t* desc, md_gto_iso_gpu_mesh_t* out) {
    const double h = desc ? iso_lattice_float(desc->cell_size) : 0.0;
    if (!s || !set || !desc || !out || !(h > 0.0)) {
        MD_LOG_ERROR("md_gto_iso_gpu_extract: invalid arguments");
        return false;
    }
    if (!isog_kernels_ready()) return false;
    {   // EXPERIMENT: the cutoff and the screening as switched
        const int screen = isog_coef_screen && set->orbitals_only && set->sh_cmax;
        if (set->eps_exp != isog_eval_eps_exp || set->screen != screen) isog_set_eps(s, (md_gto_iso_gpu_set_t*)set, isog_eval_eps_exp, screen);
    }
    memset(out, 0, sizeof(*out));
    out->num_fields = set->num_fields;
    isog_prof_begin(s);
    const md_tick_t t_start = md_tick_now();
    md_tick_t ts = t_start;
    const uint32_t F = set->num_fields;

    // root: the union of every field's surface box, as on the CPU (the same lattice for one field)
    double lo[3] = {DBL_MAX, DBL_MAX, DBL_MAX}, hi[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
    isog_field_t fields[MD_GTO_ISO_GPU_MAX_FIELDS];
    memset(fields, 0, sizeof(fields));     // all 16 uploaded: unrolled loops may read past num_fields
    uint32_t root_mask = 0;
    for (uint32_t k = 0; k < F; ++k) {
        if (!(desc->iso[k] > 0.0)) { MD_LOG_ERROR("md_gto_iso_gpu_extract: iso[%u] must be positive", k); return false; }
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
    if (L > 8) { MD_LOG_ERROR("md_gto_iso_gpu_extract: lattice too large for cell size %g", h); return false; }
    isog_mark(s, "host: root box", F, (uint32_t)L);

    memcpy(out->origin, origin, sizeof(origin));
    out->cell_size = h;
    out->cells_per_axis = 4u << L;

    md_gpu_temp_t scope = md_gpu_temp_begin(s);
    bool ok = true;
    isog_args_t a = {0};
    a.origin_h   = (md_gpu_float4){(float)origin[0], (float)origin[1], (float)origin[2], (float)h};
    a.num_shells = set->num_shells;
    a.num_ao     = set->num_ao;
    a.num_fn     = set->num_fn;
    a.fn_stride  = set->fn_stride;
    a.num_fields = F;
    a.brick_level = (uint32_t)L;
    a.step_max   = (float)(0.5 * h);
    a.lambda     = 0.5f;
    a.tol        = (float)(pow(10.0, -(double)isog_newton_tol_exp) * h);
    a.shells = set->shells.gpu; a.prims = set->prims.gpu; a.crit = set->crit.gpu;
    a.ao_nrm = set->ao_nrm.gpu; a.ao_ijk = set->ao_ijk.gpu; a.C = set->C.gpu; a.terms = set->terms.gpu;
    a.fn_tail = set->fn_tail.gpu;
    a.wts = set->wts.gpu;
    a.canon = (set->canon && isog_canon_terms) ? 1u : 0u;
    a.dmat = set->dmat.gpu;
    a.fdmat = set->fdmat.gpu;
    a.num_dmat = set->num_dmat;
    a.dform = ((isog_dform || set->factorless) && set->num_dmat > 0) ? 1u : 0u;
    // the corners in the D form: with the D-form evaluation, or alone (GTO_ISO_DFORM_CORNERS)
    const bool cdform = a.dform || (isog_dform_corners && set->num_dmat > 0 && !set->orbitals_only);
    a.diag = (uint32_t)isog_diag;
    a.drs = set->drs.gpu;
    a.dprop = set->dprop.gpu;
    // the D form's octree bounds the matrix's density, within the factors' slack of what the factored
    // evaluation gives (the slack is in the bound): either evaluation goes with it
    const bool odform = set->num_dmat > 0 && (isog_dform_octree || set->factorless);
    a.fields = isog_temp(s, sizeof(fields));
    md_gpu_upload(s, a.fields, fields, sizeof(fields));
    a.counters = isog_temp(s, sizeof(uint32_t) * 8);
    a.grid_out = isog_temp(s, sizeof(uint32_t) * 4);
    md_gpu_memset(s, a.counters, 0, sizeof(uint32_t) * 8);

    // ---- octree, one launch per level sized by the previous level's count; retried with more room
    uint32_t cap = 1u << 18, nb = 0, ctr[8] = {0};
    md_gpu_addr_t bricks = 0;
    for (int attempt = 0; attempt < 4; ++attempt) {
        md_gpu_addr_t nodes[2] = {isog_temp(s, sizeof(md_gpu_uint4) * cap), isog_temp(s, sizeof(md_gpu_uint4) * cap)};
        md_gpu_addr_t counts = isog_temp(s, sizeof(uint32_t) * (L + 2));
        if (!nodes[0] || !nodes[1] || !counts) { ok = false; break; }
        md_gpu_memset(s, counts, 0, sizeof(uint32_t) * (L + 2));
        md_gpu_memset(s, a.counters, 0, sizeof(uint32_t) * 8);
        // EXPERIMENT: GTO_ISO_SKIP_BRICK_TEST=1 keeps every child of the last tested level as a brick
        // (with its parent's sheet bits): the corners decide, which costs less than the bound when most
        // candidates are kept anyway. A superset of the bricks, so no cell is missed.
        const int Lt = isog_skip_brick_test ? L - 1 : L;
        // the first tested level l0, its nodes every node of the levels above (siblings together, as the
        // kernels append them), with the root's sheet bits
        const int l0 = Lt > 0 ? (isog_octree_start < Lt ? isog_octree_start : Lt) : 0;
        {
            md_gpu_uint4 lv[512];
            uint32_t nl = 1;
            lv[0] = (md_gpu_uint4){0, 0, 0, root_mask};
            for (int l = 0; l < l0; ++l) {
                md_gpu_uint4 nx[512];
                for (uint32_t i = 0; i < nl; ++i)
                    for (uint32_t c = 0; c < 8; ++c)
                        nx[8 * i + c] = (md_gpu_uint4){lv[i].x * 2 + (c & 1), lv[i].y * 2 + ((c >> 1) & 1), lv[i].z * 2 + (c >> 2), root_mask};
                nl *= 8;
                memcpy(lv, nx, sizeof(md_gpu_uint4) * nl);
            }
            md_gpu_upload(s, nodes[l0 & 1], lv, sizeof(md_gpu_uint4) * nl);
            md_gpu_upload(s, counts + sizeof(uint32_t) * l0, &nl, sizeof(nl));
        }
        isog_mark(s, "octree setup", 0, 0);
        uint32_t level_step[10];
        for (int l = l0; l <= Lt; ++l) {
            isog_args_t o = a;
            o.level = (uint32_t)l;
            o.nodes_in  = nodes[l & 1];
            o.nodes_out = nodes[(l + 1) & 1];
            o.count_in  = counts + sizeof(uint32_t) * l;
            o.count_out = counts + sizeof(uint32_t) * (l + 1);
            o.cap = cap;
            o.grid_mul = 1;
            o.grid_div = 0;
            o.mode = ((isog_corners_v2 && isog_corners_list) || cdform) ? 1u : 0u;     // the bricks keep their sibling slots
            MD_GPU_LAUNCH(s, isog_k[ISOG_K_GRID], md_gpu_grid(1, 1, 1), o);
            int ko = odform ? ISOG_K_ODFORM : (isog_octree_small && set->num_fn <= 4) ? ISOG_K_OCTREE_SMALL : ISOG_K_OCTREE;
            if (ko == ISOG_K_OCTREE_SMALL && !isog_k[ko]) ko = ISOG_K_OCTREE;
            if (ko == ISOG_K_ODFORM && l > 0 && isog_octree8 && isog_k[ISOG_K_ODFORM8]) {
                ko = ISOG_K_ODFORM8;    // a group per 8 siblings (below the root they come 8 at a time)
                o.grid_div = 8;
            }
            if (!isog_have(ko)) { ok = false; break; }
            md_gpu_launch_indirect(s, isog_k[ko], a.grid_out, &o, sizeof(o));
            isog_used(ko);
            level_step[l] = isog_mark(s, "octree level", 0, 0);
        }
        if (isog_prof) {    // nodes per level: in, and kept (children appended, 8 per kept node above the bricks)
            uint32_t cn[12] = {0};
            isog_read(s, cn, counts, sizeof(uint32_t) * (L + 2));
            for (int l = l0; l <= Lt; ++l) isog_mark_set(level_step[l], cn[l], l < L ? cn[l + 1] / 8 : cn[l + 1]);
            isog_mark_skip(s);
        }
        bricks = nodes[(Lt + 1) & 1];
        isog_read(s, &nb, counts + sizeof(uint32_t) * (Lt + 1), sizeof(nb));
        isog_read(s, ctr, a.counters, sizeof(ctr));
        isog_mark(s, "readback: brick count", nb, 0);
        if (!ctr[3] && nb <= cap) break;
        cap *= 4;
        nb = 0;
    }
    if (!ok || ctr[3]) {
        md_gpu_temp_end(s, scope);
        MD_LOG_ERROR("md_gto_iso_gpu_extract: octree overflow or out of device memory");
        return false;
    }
    isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_OCTREE, &ts);

    // ---- the brick list on the host: one corner slot per (brick, field with a sheet bit)
    md_gpu_uint4* bl = (md_gpu_uint4*)malloc(sizeof(md_gpu_uint4) * (nb + 1));
    uint32_t* slot = (uint32_t*)malloc(sizeof(uint32_t) * (nb + 1));
    isog_read(s, bl, bricks, sizeof(md_gpu_uint4) * nb);
    if (isog_ordered && nb > 1 && isog_sort_bricks(bl, nb)) md_gpu_upload(s, bricks, bl, sizeof(md_gpu_uint4) * nb);
    uint32_t nslots = 0, nkept = 0;
    for (uint32_t b = 0; b < nb; ++b) {
        slot[b] = nslots;
        const uint32_t m = bl[b].w;
        nslots += popcnt32((m | (m >> 1)) & 0x55555555u);
        nkept += m != 0;    // the bricks may keep their sibling slots: rejected ones have no sheet bits
    }
    out->num_bricks = nkept;
    out->num_corner_slots = nslots;
    isog_mark(s, "readback: bricks, slots", nb, nslots);
    if (isog_prof) {    // BENCHMARK: the AOs that reach each brick (what the corner kernels list), not timed
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
    a.bricks = bricks;
    a.count = nb;
    a.brick_slot = isog_temp(s, sizeof(uint32_t) * (nb + 1));
    md_gpu_upload(s, a.brick_slot, slot, sizeof(uint32_t) * nb);
    uint32_t hcap = 64;
    while (hcap < 2 * nb) hcap <<= 1;
    a.hash_mask = hcap - 1;
    a.hash_keys = isog_temp(s, sizeof(uint32_t) * hcap);
    a.hash_vals = isog_temp(s, sizeof(uint32_t) * hcap);
    a.corners   = isog_temp(s, sizeof(float) * 125 * ((size_t)nslots + 1));
    a.cell_vert = isog_temp(s, sizeof(uint32_t) * 128 * ((size_t)nslots + 1));
    const bool clist = (isog_corners_v2 && isog_corners_list) || cdform;    // the D form's corners are in list form only
    a.plist = clist ? isog_temp(s, sizeof(uint32_t) * 125 * ((size_t)nb + 1)) : 0;
    // per brick: counts, then offsets (ordered lists). Allocated either way: llvmpipe may load it in a
    // branch not taken
    a.scan  = isog_temp(s, sizeof(uint32_t) * ((size_t)nb + 1));
    a.bstate = isog_temp(s, sizeof(uint32_t) * ((size_t)nb + 1));     // flood: per brick, its round
    if (!a.scan || !a.bstate) ok = false;
    const bool ordered = isog_ordered != 0;
    if (!a.brick_slot || !a.hash_keys || !a.hash_vals || !a.corners || !a.cell_vert || (clist && !a.plist)) ok = false;
    if (ok && nb > 0) {
        md_gpu_memset(s, a.hash_keys, 0, sizeof(uint32_t) * hcap);
        MD_GPU_LAUNCH(s, isog_k[ISOG_K_HASH], md_gpu_grid_for(isog_k[ISOG_K_HASH], nb, 1, 1), a);
        isog_mark(s, "hash", nb, hcap);
        if (clist) {
            // the points to evaluate, listed brick by brick; then 256 of them per group
            const uint32_t nfn = set->num_fn;
            // APPROXIMATE: seed and flood (GTO_ISO_FLOOD): pass 0 the bricks' corners, then the seeds, then
            // round by round the bricks across faces with a sign change; else one pass over every brick
            bool flood = isog_flood == 2;
            if (isog_flood == 1) {
                flood = true;
                for (uint32_t k = 0; k < F && flood; ++k) {
                    const iso_density_t* dk = set->src[k] ? set->src[k]->f.dens : NULL;
                    if (!dk) { flood = false; break; }
                    for (int t = 0; t < dk->num_terms && flood; ++t) flood = dk->terms[t].a == dk->terms[t].b && dk->terms[t].w >= 0.0;
                }
            }
            uint32_t flood_rounds = 0, flood_bricks = 0, points_all = 0;
            for (uint32_t pass = 0; ok; ++pass) {
                md_gpu_memset(s, a.counters + sizeof(uint32_t) * 5, 0, sizeof(uint32_t));
                isog_args_t ap = a;
                ap.mode = isog_shared_corners ? 1u : 0u;
                if (flood) { ap.mode |= pass == 0 ? 16u : 8u; ap.level = pass; }
                if (ordered) {
                    // in the order of the bricks: the counts, their prefix sum, then the points at the offsets
                    ap.mode |= 2u;
                    MD_GPU_LAUNCH(s, isog_k[ISOG_K_POINTS], isog_groups(nb), ap);
                    MD_GPU_LAUNCH(s, isog_k[ISOG_K_SCAN], md_gpu_grid(1, 1, 1), ap);
                    isog_used(ISOG_K_SCAN);
                    ap.mode = (ap.mode & ~2u) | 4u;
                }
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_POINTS], isog_groups(nb), ap);
                isog_mark(s, "corners: point list", nb, 0);
                // one function: a point per thread (iso_corners_list_f1); more: the tiled GEMM, 128 points
                // per group (iso_gcorners_f*), or the register form, 256 points per group
                // densities as matrices: iso_dcorners, 128 points per group
                const bool one = nfn == 1 && F == 1;
                const bool gemm = (isog_gemm || cdform) && !one;
                // canonical sets of more than 64 functions: the wide GEMM, 64 points per group
                const int kw = nfn <= 128 ? ISOG_K_GCORN_W128 : nfn <= 192 ? ISOG_K_GCORN_W192 : ISOG_K_GCORN_W256;
                const bool wide = gemm && !cdform && a.canon && isog_gemm_wide && nfn > (isog_gemm_wide == 2 ? 16u : 64u) && isog_k[kw];
                static const int kdc[3] = {ISOG_K_DCORN, ISOG_K_DCORN2, ISOG_K_DCORN3};
                int dk = isog_dform_kernel;
                while (dk > 0 && !isog_k[kdc[dk]]) dk--;        // a D form kernel this device could make
                const int kc = one ? ((isog_corners_v4 && isog_k[ISOG_K_CLIST_F1V4]) ? ISOG_K_CLIST_F1V4 : ISOG_K_CLIST_F1)
                             : cdform ? kdc[dk]
                             : wide ? isog_wide_form(kw)
                             : gemm ? (nfn <= 16 ? ISOG_K_GCORN_F16 : (nfn <= 32 || isog_gemm_t32) ? ISOG_K_GCORN_F32
                                       : isog_gemm_pipe == 1 ? ISOG_K_GCORN_F64_P1 : isog_gemm_pipe == 2 ? ISOG_K_GCORN_F64_P2 : ISOG_K_GCORN_F64)
                             : (nfn <= 8 ? ISOG_K_CLIST_F8 : nfn <= 16 ? ISOG_K_CLIST_F16 : nfn <= 32 ? ISOG_K_CLIST_F32 : ISOG_K_CLIST_F64);
                if (!isog_have(kc)) ok = false;
                isog_args_t ac = a;
                ac.count_in = a.counters + sizeof(uint32_t) * 5;
                ac.cap = 0;
                ac.grid_mul = 1;
                ac.grid_div = wide ? 64 : gemm ? 128 : 256;
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_GRID], md_gpu_grid(1, 1, 1), ac);
                md_gpu_launch_indirect(s, isog_k[kc], a.grid_out, &ac, sizeof(ac));
                isog_used(kc);
                uint32_t fc = 0;
                const char* kn = "";
                switch (kc) {
                    case ISOG_K_CLIST_F1:  fc = 1;  kn = "corners (list f1)"; break;
                    case ISOG_K_CLIST_F1V4: fc = 1; kn = "corners (list f1 x4)"; break;
                    case ISOG_K_CLIST_F8:  fc = 8;  kn = "corners (list f8)"; break;
                    case ISOG_K_CLIST_F16: fc = 16; kn = "corners (list f16)"; break;
                    case ISOG_K_CLIST_F32: fc = 32; kn = "corners (list f32)"; break;
                    case ISOG_K_CLIST_F64: fc = 64; kn = "corners (list f64)"; break;
                    case ISOG_K_GCORN_F16: fc = 16; kn = "corners (gemm f16)"; break;
                    case ISOG_K_GCORN_F32: fc = 32; kn = "corners (gemm f32)"; break;
                    case ISOG_K_DCORN:     fc = 0;  kn = "corners (dform)"; break;
                    case ISOG_K_DCORN2:    fc = 0;  kn = "corners (dform2)"; break;
                    case ISOG_K_DCORN3:    fc = 0;  kn = "corners (dform3)"; break;
                    case ISOG_K_GCORN_F64_P1: fc = 64; kn = "corners (gemm f64 p1)"; break;
                    case ISOG_K_GCORN_F64_P2: fc = 64; kn = "corners (gemm f64 p2)"; break;
                    case ISOG_K_GCORN_W128: fc = 128; kn = "corners (gemm w128)"; break;
                    case ISOG_K_GCORN_W192: fc = 192; kn = "corners (gemm w192)"; break;
                    case ISOG_K_GCORN_W256: fc = 256; kn = "corners (gemm w256)"; break;
                    case ISOG_K_GCORN_X128: fc = 128; kn = "corners (gemm x128)"; break;
                    case ISOG_K_GCORN_X192: fc = 192; kn = "corners (gemm x192)"; break;
                    case ISOG_K_GCORN_X256: fc = 256; kn = "corners (gemm x256)"; break;
                    case ISOG_K_GCORN_Y128: fc = 128; kn = "corners (gemm y128)"; break;
                    case ISOG_K_GCORN_Y192: fc = 192; kn = "corners (gemm y192)"; break;
                    case ISOG_K_GCORN_Y256: fc = 256; kn = "corners (gemm y256)"; break;
                    default:               fc = 64; kn = "corners (gemm f64)"; break;
                }
                if (isog_prof) isog_prof->corner_cols = fc ? fc * ((nfn + fc - 1) / fc) : 0;
                const uint32_t ck = isog_mark(s, kn, nkept, nslots);
                if (isog_prof && flood) {   // aux: the lattice points evaluated, summed over the passes
                    uint32_t c8[8];
                    isog_read(s, c8, a.counters, sizeof(c8));
                    points_all += c8[5];
                    isog_mark_set(ck, nkept, c8[5]);
                    isog_prof->corner_points = points_all;
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
                        isog_prof->cgroup_size = wide ? 64u : gemm ? 128u : 256u;
                        isog_prof_groups(set, px, 3, np, isog_prof->cgroup_size, &isog_prof->cgroup_ao_sum, &isog_prof->cgroup_count);
                    }
                    free(pl);
                    free(px);
                    isog_mark_skip(s);
                }
                if (isog_diag && !flood) md_gpu_memset(s, a.corners, 0, sizeof(float) * 125 * (size_t)nslots);   // DIAGNOSTIC: no surface
                if (isog_shared_corners) {
                    isog_args_t af = a;
                    af.mode = flood ? (pass == 0 ? 16u : 8u) : 0u;
                    af.level = pass;
                    MD_GPU_LAUNCH(s, isog_k[ISOG_K_CORNERS_FILL], isog_groups(nb), af);
                    isog_mark(s, "corners: fill shared", nb, 0);
                }
                if (!flood) break;
                // the next round: the seeds (pass 0), or the bricks across this round's faces with a sign change
                md_gpu_memset(s, a.counters + sizeof(uint32_t) * 6, 0, sizeof(uint32_t));
                isog_args_t ar = a;
                ar.level = pass;
                if (pass == 0) {
                    MD_GPU_LAUNCH(s, isog_k[ISOG_K_FSEED], md_gpu_grid((nb + 63) / 64, 1, 1), ar);
                    isog_used(ISOG_K_FSEED);
                } else {
                    MD_GPU_LAUNCH(s, isog_k[ISOG_K_FFACE], isog_groups(nb), ar);
                    isog_used(ISOG_K_FFACE);
                }
                uint32_t next = 0;
                isog_read(s, &next, a.counters + sizeof(uint32_t) * 6, sizeof(next));
                isog_mark(s, pass == 0 ? "flood: seeds" : "flood: faces", nb, next);
                flood_bricks += next;
                if (next == 0 || pass >= 64) break;
                flood_rounds++;
            }
            if (flood && ok) {
                // the bricks no round reached: no cells (and the host's list says so, for the profile)
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_FFINISH], md_gpu_grid((nb + 63) / 64, 1, 1), a);
                isog_used(ISOG_K_FFINISH);
                isog_read(s, bl, a.bricks, sizeof(md_gpu_uint4) * nb);
                isog_mark(s, "flood: finish", flood_bricks, flood_rounds);
                if (isog_diag) md_gpu_memset(s, a.corners, 0, sizeof(float) * 125 * (size_t)nslots);   // DIAGNOSTIC: no surface
            }
        } else if (isog_corners_v2) {
            // the function columns in registers: all of the set's, or tiles of 64
            const uint32_t nfn = set->num_fn;
            const int kc = (nfn == 1 && F == 1) ? ISOG_K_CORNERS_F1 : nfn <= 8 ? ISOG_K_CORNERS_F8 : nfn <= 16 ? ISOG_K_CORNERS_F16
                         : nfn <= 32 ? ISOG_K_CORNERS_F32 : ISOG_K_CORNERS_F64;
            isog_args_t ac = a;
            ac.mode = isog_shared_corners ? 1u : 0u;
            MD_GPU_LAUNCH(s, isog_k[kc], isog_groups(nb), ac);
            isog_used(kc);
            static const uint32_t kcols[] = {1, 8, 16, 32, 64};
            if (isog_prof) isog_prof->corner_cols = kcols[kc - ISOG_K_CORNERS_F1] * ((nfn + kcols[kc - ISOG_K_CORNERS_F1] - 1) / kcols[kc - ISOG_K_CORNERS_F1]);
            static const char* kname[] = {"corners (f1)", "corners (f8)", "corners (f16)", "corners (f32)", "corners (f64)"};
            const uint32_t ck = isog_mark(s, kname[kc - ISOG_K_CORNERS_F1], nkept, nslots);
            if (isog_prof) {    // aux: the lattice points evaluated (of 125 per brick)
                uint32_t c8[8];
                isog_read(s, c8, a.counters, sizeof(c8));
                isog_mark_set(ck, nkept, c8[5]);
                isog_prof->corner_points = c8[5];
                isog_mark_skip(s);
            }
            if (isog_shared_corners) {
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_CORNERS_FILL], isog_groups(nb), a);
                isog_mark(s, "corners: fill shared", nb, 0);
            }
        } else {
            MD_GPU_LAUNCH(s, isog_k[set->orbitals_only ? ISOG_K_CORNERS1 : ISOG_K_CORNERS], isog_groups(nb), a);
            isog_used(set->orbitals_only ? ISOG_K_CORNERS1 : ISOG_K_CORNERS);
            if (isog_prof) isog_prof->corner_cols = set->orbitals_only ? set->num_fn : set->fn_stride;
            isog_mark(s, set->orbitals_only ? "corners (orbital)" : "corners (gemm)", nb, nslots);
        }
    }
    isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_CORNERS, &ts);

    // ---- Surface Nets: count, then size the vertex buffers and emit
    uint32_t nv = 0, nq = 0;
    if (ok && nb > 0) {
        a.mode = ordered ? 2u : 0u;     // bit 2: the brick's count into scan, for vertices in the order of the bricks
        md_gpu_memset(s, a.counters, 0, sizeof(uint32_t) * 7);     // [7] (dform overflow) kept
        MD_GPU_LAUNCH(s, isog_k[ISOG_K_CELLS], isog_groups(nb), a);
        if (ordered) {
            MD_GPU_LAUNCH(s, isog_k[ISOG_K_SCAN], md_gpu_grid(1, 1, 1), a);
            isog_used(ISOG_K_SCAN);
        }
        isog_mark(s, "cells: count", nb, 0);
        isog_read(s, &nv, a.counters, sizeof(nv));
        isog_mark(s, "readback: vertex count", nv, 0);
        if (isog_prof && ordered) {     // BENCHMARK: the bricks given corners, by their corners (not timed)
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
        memcpy(dbg->origin, origin, sizeof(origin));
        dbg->brick_level = L;
        isog_mark_skip(s);
    }
    free(bl);
    free(slot);

    const uint32_t cap_v = nv + 1, cap_q = 3 * nv + 1;
    a.vpos     = isog_temp(s, sizeof(md_gpu_float4) * cap_v);
    a.vpos_out = isog_temp(s, sizeof(md_gpu_float4) * cap_v);
    a.vnrm     = isog_temp(s, sizeof(md_gpu_float4) * cap_v);
    a.vval     = isog_temp(s, sizeof(md_gpu_float4) * cap_v);
    a.vinfo    = isog_temp(s, sizeof(md_gpu_uint4) * cap_v);
    a.quads    = isog_temp(s, sizeof(md_gpu_uint4) * cap_q);
    a.remap    = isog_temp(s, sizeof(uint32_t) * cap_v);
    a.bucket   = isog_temp(s, sizeof(uint32_t) * 192);
    a.vdone    = isog_temp(s, sizeof(uint32_t) * cap_v);
    if (!a.vpos || !a.vpos_out || !a.vnrm || !a.vval || !a.vinfo || !a.quads || !a.remap || !a.bucket || !a.vdone) ok = false;
    if (ok && nv > 0) {
        md_gpu_memset(s, a.cell_vert, 0xFF, sizeof(uint32_t) * 128 * (size_t)nslots);
        md_gpu_memset(s, a.counters, 0, sizeof(uint32_t) * 7);     // [7] (dform overflow) kept
        a.mode = ordered ? 3u : 1u;
        a.cap = cap_v;
        MD_GPU_LAUNCH(s, isog_k[ISOG_K_CELLS], isog_groups(nb), a);
        isog_mark(s, "cells: emit", nb, nv);
        a.cap = cap_q;
        MD_GPU_LAUNCH(s, isog_k[ISOG_K_QUADS], isog_groups(nb), a);
        isog_mark(s, "quads", nb, 0);
    }
    isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_NETS, &ts);

    // ---- vertices: Newton projection, relaxation; the quad count is read with the final counts
    if (ok && nv > 0) {
        isog_args_t v = a;
        v.count = nv;
        const md_gpu_grid_t g_vert = md_gpu_grid_for(isog_k[ISOG_K_NEWTON], nv, 1, 1);
        int ke = set->orbitals_only ? ISOG_K_EVAL1 : ISOG_K_EVAL;
        const bool geval = !set->orbitals_only && ((isog_eval_v2 && isog_gemm) || a.dform);
        const uint32_t nfv = set->num_fn;
        const int kwv = nfv <= 128 ? ISOG_K_GEVAL_W128 : nfv <= 192 ? ISOG_K_GEVAL_W192 : ISOG_K_GEVAL_W256;
        if (a.dform) {
            static const int kde[3] = {ISOG_K_DEVAL, ISOG_K_DEVAL2, ISOG_K_DEVAL3};
            int dk = isog_dform_kernel;
            while (dk > 0 && !isog_k[kde[dk]]) dk--;
            ke = kde[dk];
        } else if (geval && a.canon && isog_gemm_wide && nfv > (isog_gemm_wide == 2 ? 16u : 64u) && isog_k[kwv]) {
            ke = isog_wide_form(kwv);
        } else if (geval) {
            const uint32_t nfn = set->num_fn;
            ke = nfn <= 16 ? ISOG_K_GEVAL_F16 : (nfn <= 32 || isog_gemm_t32) ? ISOG_K_GEVAL_F32
               : isog_gemm_pipe == 1 ? ISOG_K_GEVAL_F64_P1 : isog_gemm_pipe == 2 ? ISOG_K_GEVAL_F64_P2 : ISOG_K_GEVAL_F64;
        } else if (!set->orbitals_only && isog_eval_v2) {
            const uint32_t nfn = set->num_fn;
            ke = nfn <= 8 ? ISOG_K_EVAL_F8 : nfn <= 16 ? ISOG_K_EVAL_F16 : nfn <= 32 ? ISOG_K_EVAL_F32 : ISOG_K_EVAL_F64;
        }
        const bool elist = set->orbitals_only && set->num_fn == 1 && set->num_fields == 1 && isog_eval_list;
        if (elist) ke = ISOG_K_EVAL_LIST1;
        const bool vwide = (ke >= ISOG_K_GEVAL_W128 && ke <= ISOG_K_GEVAL_W256) || (ke >= ISOG_K_GEVAL_X128 && ke <= ISOG_K_GEVAL_X256) ||
                           (ke >= ISOG_K_GEVAL_Y128 && ke <= ISOG_K_GEVAL_Y256);
        if (!isog_have(ke)) ok = false;
        md_gpu_kernel_t k_eval = isog_k[ke];
        isog_used(ke);
        if (isog_prof) {    // BENCHMARK: the boxes of the vertex groups (consecutive vertices), not timed
            const uint32_t vg = elist ? 128u : vwide ? 16u : geval ? 32u : set->orbitals_only ? 0u : 64u;
            md_gpu_float4* vp = vg ? (md_gpu_float4*)malloc(sizeof(md_gpu_float4) * nv) : NULL;
            if (vp && isog_read(s, vp, a.vpos, sizeof(md_gpu_float4) * nv)) {
                isog_prof->vgroup_size = vg;
                isog_prof_groups(set, &vp[0].x, 4, nv, vg, &isog_prof->vgroup_ao_sum, &isog_prof->vgroup_count);
            }
            free(vp);
            isog_mark_skip(s);
        }
        // groups of 32 vertices (the GEMM: 4 channels each, 128 rows), 64 (the register form), 128 (one
        // orbital from lists) or a thread per vertex
        const md_gpu_grid_t g_eval = elist ? isog_groups((nv + 127) / 128) : set->orbitals_only ? g_vert
                                   : isog_groups(vwide ? (nv + 15) / 16 : geval ? (nv + 31) / 32 : (nv + 63) / 64);
        // Newton steps; converged vertices are flagged and skipped (they no longer move, so their
        // values stay valid for the final normals)
        #define ISOG_PROJECT(iters) do {                                               \
            md_gpu_memset(s, v.vdone, 0, sizeof(uint32_t) * nv);                       \
            for (uint32_t it = 0; it <= (iters); ++it) {                               \
                uint32_t act_, grp_;                                                   \
                isog_prof_active(s, v.vdone, nv, &act_, &grp_);                        \
                MD_GPU_LAUNCH(s, k_eval, g_eval, v);                                   \
                isog_mark(s, "vertex eval", act_, grp_);                               \
                v.mode = it == (iters) ? 1u : (isog_newton_predict ? 2u : 0u);         \
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_NEWTON], g_vert, v);                    \
                isog_mark(s, v.mode == 1u ? "vertex normals" : "vertex newton", act_, 0);\
            }                                                                          \
        } while (0)
        if (isog_relax_first) {
            // relaxed on the lattice's own (trilinear) surface, with its normals from iso_cells, then
            // projected once: no field evaluation before the projection
            for (uint32_t r = 0; r < desc->relax_iterations; ++r) {
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_RELAX], g_vert, v);
                isog_mark(s, "relax", nv, 0);
                md_gpu_copy(s, v.vpos, v.vpos_out, sizeof(md_gpu_float4) * nv);
                isog_mark(s, "relax: copy", nv, 0);
            }
            isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_RELAX, &ts);
            ISOG_PROJECT(desc->project_iterations);
            for (int r = 1; r < isog_relax_first && desc->relax_iterations; ++r) {   // EXPERIMENT: polish with the true normals
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_RELAX], g_vert, v);
                isog_mark(s, "relax", nv, 0);
                md_gpu_copy(s, v.vpos, v.vpos_out, sizeof(md_gpu_float4) * nv);
                isog_mark(s, "relax: copy", nv, 0);
                ISOG_PROJECT(3u);
            }
            isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_PROJECT, &ts);
        } else {
            ISOG_PROJECT(desc->project_iterations);
            isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_PROJECT, &ts);
            for (uint32_t r = 0; r < desc->relax_iterations; ++r) {
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_RELAX], g_vert, v);
                isog_mark(s, "relax", nv, 0);
                md_gpu_copy(s, v.vpos, v.vpos_out, sizeof(md_gpu_float4) * nv);
                isog_mark(s, "relax: copy", nv, 0);
                ISOG_PROJECT(3u);
            }
            isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_RELAX, &ts);
        }
        #undef ISOG_PROJECT
    }

    // ---- regroup per (field, sheet) into the output buffers
    if (ok && nv > 0) {
        isog_read(s, ctr, a.counters, sizeof(ctr));
        nq = ctr[1] < cap_q ? ctr[1] : cap_q;
        out->open_edges = ctr[2];
        if (ctr[3]) ok = false;
        if (ctr[7]) {
            MD_LOG_ERROR("md_gto_iso_gpu_extract: more than 2048 shells (or 8192 AOs) reach a group of points (the D form)");
            ok = false;
        }
        isog_mark(s, "readback: counters", nv, nq);
    }
    if (ok) {
        md_gpu_mem_t pos = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(md_gpu_float4) * cap_v);
        md_gpu_mem_t nrm = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(md_gpu_float4) * cap_v);
        md_gpu_mem_t idx = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(uint32_t) * 6 * ((size_t)nq + 1));
        out->position = pos.gpu; out->normal = nrm.gpu; out->index = idx.gpu;
        if (!pos.gpu || !nrm.gpu || !idx.gpu) ok = false;
    }
    if (ok && nv > 0) {
        a.out_pos = out->position; a.out_nrm = out->normal; a.out_idx = out->index;
        a.count = nv;
        a.cap = nq;
        isog_mark(s, "output alloc", nv, nq);
        md_gpu_memset(s, a.bucket, 0, sizeof(uint32_t) * 192);
        MD_GPU_LAUNCH(s, isog_k[ISOG_K_BCOUNT], md_gpu_grid_for(isog_k[ISOG_K_BCOUNT], nv > nq ? nv : nq, 1, 1), a);
        MD_GPU_LAUNCH(s, isog_k[ISOG_K_BSCAN], md_gpu_grid(1, 1, 1), a);
        isog_mark(s, "bucket: count, scan", nv, nq);
        MD_GPU_LAUNCH(s, isog_k[ISOG_K_BVERTS], md_gpu_grid_for(isog_k[ISOG_K_BVERTS], nv, 1, 1), a);
        isog_mark(s, "bucket: vertices", nv, 0);
        if (nq) MD_GPU_LAUNCH(s, isog_k[ISOG_K_BTRIS], md_gpu_grid_for(isog_k[ISOG_K_BTRIS], nq, 1, 1), a);
        isog_mark(s, "bucket: triangles", nq, 0);
        uint32_t bk[128];
        isog_read(s, bk, a.bucket, sizeof(bk));
        isog_mark(s, "readback: buckets", 0, 0);
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
        out->index_count = 6 * nq;
    }
    isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_OUTPUT, &ts);
    md_gpu_temp_end(s, scope);
    out->seconds_total = md_tick_to_seconds(md_tick_now() - t_start);
    if (!ok) {
        md_gto_iso_gpu_mesh_free(s, out);   // keeps nothing, not even the lattice
        MD_LOG_ERROR("md_gto_iso_gpu_extract: out of device memory or buffer overflow");
    }
    return ok;
}

void md_gto_iso_gpu_mesh_free(md_gpu_stream_t stream, md_gto_iso_gpu_mesh_t* mesh) {
    if (!mesh) return;
    md_gpu_free(stream, mesh->position);
    md_gpu_free(stream, mesh->normal);
    md_gpu_free(stream, mesh->index);
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

// md_gto_iso_gpu.inl - the md_gpu path of md_gto_iso (included by md_gto_iso.c under MD_ENABLE_GPU).
//
// Host side of src/shaders/iso/gto_iso_gpu.slang: builds a field set (basis, coefficient matrix C,
// terms per field), then issues the pipeline: octree level by level, bricks, corner values, Surface
// Nets, Newton projection and relaxation, per-field regrouping. Three waits for the GPU: the brick
// list (to size the corner storage), the vertex count (to size the mesh buffers), the final counts.

#include <core/md_intrinsics.h>
#include <core/md_gpu.h>
#include <gto_iso_gpu_shaders.inl>

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
    md_gpu_addr_t shells, prims, crit, ao_nrm, ao_ijk, C, terms, fields;
    md_gpu_addr_t nodes_in, nodes_out, count_in, count_out, grid_out;
    md_gpu_addr_t bricks, brick_slot, hash_keys, hash_vals, corners, cell_vert;
    md_gpu_addr_t vpos, vpos_out, vnrm, vval, vinfo, quads, counters, bucket, remap, vdone;
    md_gpu_addr_t out_pos, out_nrm, out_idx;
    md_gpu_addr_t fn_tail;
} isog_args_t;

enum {
    ISOG_K_GRID, ISOG_K_OCTREE, ISOG_K_HASH, ISOG_K_CORNERS, ISOG_K_EVAL, ISOG_K_CELLS, ISOG_K_QUADS,
    ISOG_K_NEWTON, ISOG_K_RELAX, ISOG_K_BCOUNT, ISOG_K_BSCAN, ISOG_K_BVERTS, ISOG_K_BTRIS, ISOG_K_EVAL1, ISOG_K_CORNERS1,
    ISOG_K_CORNERS_F1, ISOG_K_CORNERS_F8, ISOG_K_CORNERS_F16, ISOG_K_CORNERS_F32, ISOG_K_CORNERS_F64,
    ISOG_K_EVAL_F8, ISOG_K_EVAL_F16, ISOG_K_EVAL_F32, ISOG_K_EVAL_F64, ISOG_K_CORNERS_FILL, ISOG_K_COUNT
};

static md_gpu_kernel_t isog_k[ISOG_K_COUNT];

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
}

static bool isog_kernels_ready(void) {
    if (isog_k[ISOG_K_COUNT - 1]) return true;
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
    };
    for (int i = 0; i < ISOG_K_COUNT; ++i) {
        if (!isog_k[i]) isog_k[i] = md_gpu_kernel_create(isog_device, &d[i]);
        if (!isog_k[i]) {
            MD_LOG_ERROR("md_gto_iso_gpu: kernel %d: %s", i, md_gpu_last_error() ? md_gpu_last_error() : "?");
            isog_failed = true;
            return false;
        }
    }
    return true;
}

void md_gto_iso_gpu_shutdown(void) {
    for (int i = 0; i < ISOG_K_COUNT; ++i) {
        if (isog_k[i]) { md_gpu_kernel_destroy(isog_k[i]); isog_k[i] = NULL; }
    }
    isog_device = NULL;
    isog_failed = false;
}

// ---------------------------------------------------------------------------------------------
// The field set
// ---------------------------------------------------------------------------------------------

#define ISOG_FN_CAP 1024     // O_FN_CAP in the shader

struct md_gto_iso_gpu_set_t {
    md_gpu_mem_t shells, prims, crit, ao_nrm, ao_ijk, C, terms, fn_tail;
    uint32_t     num_shells, num_ao, num_fn, fn_stride, num_fields;
    bool         orbitals_only;     // every field one function, one linear term: iso_eval_orbital
    isog_field_t fields[MD_GTO_ISO_GPU_MAX_FIELDS];    // iso filled in per extraction
    iso_field_t  host[MD_GTO_ISO_GPU_MAX_FIELDS];         // shallow copies for the root box (their tails): the fields outlive the set
    const md_gto_iso_field_t* src[MD_GTO_ISO_GPU_MAX_FIELDS];
    struct md_allocator_i* alloc;
    isog_shell_t* prof_shells;      // BENCHMARK: the shells on the host, for the profile's work counts
};

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
    if (!ok || nf > ISOG_FN_CAP) {
        md_array_free(B.radial.shells, md_get_heap_allocator()); md_array_free(B.radial.crit, md_get_heap_allocator()); free(B.ao_shells);
        MD_LOG_ERROR("md_gto_iso_gpu_set_create: fields of a different basis, or more than %d functions in all", ISOG_FN_CAP);
        return NULL;
    }
    const uint32_t stride = ((nf + 63u) / 64u) * 64u;
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
        if (fk->f.dens) {
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

    // shells, primitives, radial critical intervals (outward rounded), AO tables
    isog_shell_t* sh = (isog_shell_t*)calloc(B.num_ao_shells + 1, sizeof(isog_shell_t));
    md_gpu_float2* pr = NULL;
    md_gpu_float4* cr = NULL;
    float*    nrm = (float*)calloc(nao + 1, sizeof(float));
    uint32_t* ijk = (uint32_t*)calloc(nao + 1, sizeof(uint32_t));
    for (int s = 0; s < B.num_ao_shells; ++s) {
        const iso_ao_shell_t* a = &B.ao_shells[s];
        const iso_shell_t* g = &B.radial.shells[s];
        // point values: every AO of the shell below 1e-10 beyond this radius (the corner values only
        // decide on which side of an iso value a corner lies; the bounds do not use it)
        double nmax = 0.0, amin = DBL_MAX;
        for (int n = 0; n < a->ncart; ++n) nmax = fmax(nmax, fabs(a->nrm[n]));
        for (int q = 0; q < a->nprim; ++q) amin = fmin(amin, a->alpha[q]);
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
    set->num_shells = (uint32_t)B.num_ao_shells;
    set->num_ao = nao;
    set->num_fn = nf;
    set->fn_stride = stride;
    set->num_fields = num_fields;
    set->orbitals_only = orbitals_only;

    // per function: what the shells that do not reach a box can add over it, every AO of theirs being
    // below 1e-10 there (shells' w); twice that, for the float distances that decide it
    float* tail = (float*)calloc(stride + 1, sizeof(float));
    for (uint32_t j = 0; j < nf; ++j) {
        double t = 0.0;
        for (uint32_t i = 0; i < nao; ++i) t += fabs((double)C[(size_t)i * stride + j]);
        tail[j] = isog_f32_up(2.0e-10 * t * (1.0 + 1.0e-6));
    }

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
    }
    free(tail);
#undef ISOG_UP
    set->prof_shells = sh;          // BENCHMARK: kept for the profile's work counts
    free(nrm); free(ijk); free(C);
    md_array_free(pr, heap); md_array_free(cr, heap); md_array_free(terms, heap);
    md_array_free(B.radial.shells, heap); md_array_free(B.radial.crit, heap); free(B.ao_shells);
    if (!ok) {
        MD_LOG_ERROR("md_gto_iso_gpu_set_create: a pair term spans two tiles of 64 functions, or out of device memory");
        md_gto_iso_gpu_set_destroy(stream, set);
        return NULL;
    }
    return set;
}

void md_gto_iso_gpu_set_destroy(md_gpu_stream_t stream, md_gto_iso_gpu_set_t* set) {
    if (!set) return;
    md_gpu_free(stream, set->shells.gpu); md_gpu_free(stream, set->prims.gpu); md_gpu_free(stream, set->crit.gpu);
    md_gpu_free(stream, set->ao_nrm.gpu); md_gpu_free(stream, set->ao_ijk.gpu); md_gpu_free(stream, set->C.gpu);
    md_gpu_free(stream, set->terms.gpu); md_gpu_free(stream, set->fn_tail.gpu);
    free(set->prof_shells);
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

const char* md_gto_iso_gpu_profile_variant(void) {
    static char buf[128];
    snprintf(buf, sizeof(buf), "GTO_ISO_RELAX_FIRST=%d GTO_ISO_CORNERS_V2=%d GTO_ISO_EVAL_V2=%d GTO_ISO_NEWTON_PREDICT=%d GTO_ISO_SKIP_BRICK_TEST=%d GTO_ISO_SHARED_CORNERS=%d",
             isog_relax_first, isog_corners_v2, isog_eval_v2, isog_newton_predict, isog_skip_brick_test, isog_shared_corners);
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

static void isog_prof_begin(md_gpu_stream_t s) {
    if (!isog_prof) return;
    isog_prof->num_steps = 0;
    isog_prof->brick_ao_sum = 0;
    isog_prof->brick_ao_max = 0;
    isog_prof->corner_cols = 0;
    isog_prof->corner_points = 0;
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
            const double fma = pts * (double)p->brick_ao_sum * cols;
            ISOG_PR("    corners: %.1f AOs reach a brick (max %u of %u), %.1f points evaluated per brick, %u functions: %.3g multiply-adds, %.1f G/s\n",
                    (double)p->brick_ao_sum / (double)items[k], p->brick_ao_max, p->num_ao, pts, p->num_fn, fma, fma / (w * 1e6));
        }
        if (strcmp(name[k], "octree level") == 0 && items[k] > 0) {
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
    ISOG_PR("    octree levels, nodes in -> kept:");
    for (uint32_t i = 0; i < p->num_steps; ++i) {
        if (strcmp(p->step[i].name, "octree level") != 0) continue;
        ISOG_PR(" %u->%u", p->step[i].items, p->step[i].aux);
    }
    ISOG_PR("\n");
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
    a.tol        = (float)(1.0e-4 * h);
    a.shells = set->shells.gpu; a.prims = set->prims.gpu; a.crit = set->crit.gpu;
    a.ao_nrm = set->ao_nrm.gpu; a.ao_ijk = set->ao_ijk.gpu; a.C = set->C.gpu; a.terms = set->terms.gpu;
    a.fn_tail = set->fn_tail.gpu;
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
        const md_gpu_uint4 root = {0, 0, 0, root_mask};
        const uint32_t one = 1;
        md_gpu_upload(s, nodes[0], &root, sizeof(root));
        md_gpu_upload(s, counts, &one, sizeof(one));
        isog_mark(s, "octree setup", 0, 0);
        uint32_t level_step[10];
        // EXPERIMENT: GTO_ISO_SKIP_BRICK_TEST=1 keeps every child of the last tested level as a brick
        // (with its parent's sheet bits): the corners decide, which costs less than the bound when most
        // candidates are kept anyway. A superset of the bricks, so no cell is missed.
        const int Lt = isog_skip_brick_test ? L - 1 : L;
        for (int l = 0; l <= Lt; ++l) {
            isog_args_t o = a;
            o.level = (uint32_t)l;
            o.nodes_in  = nodes[l & 1];
            o.nodes_out = nodes[(l + 1) & 1];
            o.count_in  = counts + sizeof(uint32_t) * l;
            o.count_out = counts + sizeof(uint32_t) * (l + 1);
            o.cap = cap;
            o.grid_mul = 1;
            o.grid_div = 0;
            MD_GPU_LAUNCH(s, isog_k[ISOG_K_GRID], md_gpu_grid(1, 1, 1), o);
            md_gpu_launch_indirect(s, isog_k[ISOG_K_OCTREE], a.grid_out, &o, sizeof(o));
            level_step[l] = isog_mark(s, "octree level", 0, 0);
        }
        if (isog_prof) {    // nodes per level: in, and kept (children appended, 8 per kept node above the bricks)
            uint32_t cn[12] = {0};
            isog_read(s, cn, counts, sizeof(uint32_t) * (L + 2));
            for (int l = 0; l <= Lt; ++l) isog_mark_set(level_step[l], cn[l], l < L ? cn[l + 1] / 8 : cn[l + 1]);
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
    out->num_bricks = nb;
    isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_OCTREE, &ts);

    // ---- the brick list on the host: one corner slot per (brick, field with a sheet bit)
    md_gpu_uint4* bl = (md_gpu_uint4*)malloc(sizeof(md_gpu_uint4) * (nb + 1));
    uint32_t* slot = (uint32_t*)malloc(sizeof(uint32_t) * (nb + 1));
    isog_read(s, bl, bricks, sizeof(md_gpu_uint4) * nb);
    uint32_t nslots = 0;
    for (uint32_t b = 0; b < nb; ++b) {
        slot[b] = nslots;
        const uint32_t m = bl[b].w;
        nslots += popcnt32((m | (m >> 1)) & 0x55555555u);
    }
    out->num_corner_slots = nslots;
    isog_mark(s, "readback: bricks, slots", nb, nslots);
    if (isog_prof) {    // BENCHMARK: the AOs that reach each brick (what the corner kernels list), not timed
        uint64_t sum = 0;
        uint32_t mx = 0;
        const float he = (float)h;
        for (uint32_t b = 0; b < nb; ++b) {
            const float blo[3] = {(float)origin[0] + (float)(bl[b].x * 4u) * he, (float)origin[1] + (float)(bl[b].y * 4u) * he, (float)origin[2] + (float)(bl[b].z * 4u) * he};
            uint32_t n = 0;
            for (uint32_t k = 0; k < set->num_shells; ++k) {
                const isog_shell_t* sh = &set->prof_shells[k];
                const float c[3] = {sh->c.x, sh->c.y, sh->c.z};
                float d2 = 0.0f;
                for (int ax = 0; ax < 3; ++ax) {
                    const float d = fmaxf(fmaxf(blo[ax] - c[ax], c[ax] - (blo[ax] + 4.0f * he)), 0.0f);
                    d2 += d * d;
                }
                if (d2 <= sh->c.w) n += sh->info.y;
            }
            sum += n;
            mx = n > mx ? n : mx;
        }
        isog_prof->brick_ao_sum = sum;
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
    if (!a.brick_slot || !a.hash_keys || !a.hash_vals || !a.corners || !a.cell_vert) ok = false;
    if (ok && nb > 0) {
        md_gpu_memset(s, a.hash_keys, 0, sizeof(uint32_t) * hcap);
        MD_GPU_LAUNCH(s, isog_k[ISOG_K_HASH], md_gpu_grid_for(isog_k[ISOG_K_HASH], nb, 1, 1), a);
        isog_mark(s, "hash", nb, hcap);
        if (isog_corners_v2) {
            // the function columns in registers: all of the set's, or tiles of 64
            const uint32_t nfn = set->num_fn;
            const int kc = (nfn == 1 && F == 1) ? ISOG_K_CORNERS_F1 : nfn <= 8 ? ISOG_K_CORNERS_F8 : nfn <= 16 ? ISOG_K_CORNERS_F16
                         : nfn <= 32 ? ISOG_K_CORNERS_F32 : ISOG_K_CORNERS_F64;
            isog_args_t ac = a;
            ac.mode = isog_shared_corners ? 1u : 0u;
            MD_GPU_LAUNCH(s, isog_k[kc], isog_groups(nb), ac);
            static const uint32_t kcols[] = {1, 8, 16, 32, 64};
            if (isog_prof) isog_prof->corner_cols = kcols[kc - ISOG_K_CORNERS_F1] * ((nfn + kcols[kc - ISOG_K_CORNERS_F1] - 1) / kcols[kc - ISOG_K_CORNERS_F1]);
            static const char* kname[] = {"corners (f1)", "corners (f8)", "corners (f16)", "corners (f32)", "corners (f64)"};
            const uint32_t ck = isog_mark(s, kname[kc - ISOG_K_CORNERS_F1], nb, nslots);
            if (isog_prof) {    // aux: the lattice points evaluated (of 125 per brick)
                uint32_t c8[8];
                isog_read(s, c8, a.counters, sizeof(c8));
                isog_mark_set(ck, nb, c8[5]);
                isog_prof->corner_points = c8[5];
                isog_mark_skip(s);
            }
            if (isog_shared_corners) {
                MD_GPU_LAUNCH(s, isog_k[ISOG_K_CORNERS_FILL], isog_groups(nb), a);
                isog_mark(s, "corners: fill shared", nb, 0);
            }
        } else {
            MD_GPU_LAUNCH(s, isog_k[set->orbitals_only ? ISOG_K_CORNERS1 : ISOG_K_CORNERS], isog_groups(nb), a);
            if (isog_prof) isog_prof->corner_cols = set->orbitals_only ? set->num_fn : set->fn_stride;
            isog_mark(s, set->orbitals_only ? "corners (orbital)" : "corners (gemm)", nb, nslots);
        }
    }
    isog_stage(s, desc, out, MD_GTO_ISO_GPU_STAGE_CORNERS, &ts);

    // ---- Surface Nets: count, then size the vertex buffers and emit
    uint32_t nv = 0, nq = 0;
    if (ok && nb > 0) {
        a.mode = 0;
        md_gpu_memset(s, a.counters, 0, sizeof(uint32_t) * 8);
        MD_GPU_LAUNCH(s, isog_k[ISOG_K_CELLS], isog_groups(nb), a);
        isog_mark(s, "cells: count", nb, 0);
        isog_read(s, &nv, a.counters, sizeof(nv));
        isog_mark(s, "readback: vertex count", nv, 0);
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
        md_gpu_memset(s, a.counters, 0, sizeof(uint32_t) * 8);
        a.mode = 1;
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
        md_gpu_kernel_t k_eval = set->orbitals_only ? isog_k[ISOG_K_EVAL1] : isog_k[ISOG_K_EVAL];
        if (!set->orbitals_only && isog_eval_v2) {
            const uint32_t nfn = set->num_fn;
            k_eval = isog_k[nfn <= 8 ? ISOG_K_EVAL_F8 : nfn <= 16 ? ISOG_K_EVAL_F16 : nfn <= 32 ? ISOG_K_EVAL_F32 : ISOG_K_EVAL_F64];
        }
        const md_gpu_grid_t g_eval = set->orbitals_only ? g_vert : isog_groups((nv + 63) / 64);
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

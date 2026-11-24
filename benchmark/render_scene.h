// render_scene.h: a structure through md_render_mol, as md_render_snap and
// md_render_window draw it: the system loaded and inferred, its topology and
// per-frame block written, the atoms' colours and selection uploaded, and the
// representations asked for emitted as md_render batches. Header only, for
// those two programs.

#ifndef RENDER_SCENE_H
#define RENDER_SCENE_H

#include <md_system.h>
#include <md_system_state.h>
#include <md_infer.h>
#include <md_io.h>
#include <md_element.h>
#include <md_backbone.h>
#include <md_render.h>
#include <md_render_mol.h>
#include <core/md_gpu.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct scene_rep_t {
    md_render_mol_rep_type_t type;
    float                    param;      /* < 0: the type's default */
    md_render_fade_t         fade;
    float                    residues[2];
} scene_rep_t;

typedef struct scene_desc_t {
    const char* structure;
    const char* coords;                  /* NULL: the structure's own */
    scene_rep_t reps[8];
    uint32_t    rep_count;               /* 0: spacefill */
    uint32_t    replicate;               /* N x N x N instances along the cell (0 or 1: one) */
    int64_t     sel_beg, sel_end;        /* atoms selected (the mask); sel_end <= sel_beg: none */
    bool        chain_gaps, bond_orders;
} scene_desc_t;

typedef struct scene_t {
    md_allocator_i*       alloc;
    md_system_t           sys;
    md_system_state_t     state;
    uint32_t              n;             /* atoms */
    md_render_mol_topo_t  topo;
    md_render_mol_atoms_t atoms;
    md_gpu_addr_t         dyn;
    uint64_t              content[MD_RENDER_CONTENT_COUNT];
    md_render_batch_t     batches[64];
    uint32_t              batch_count;
    uint32_t              n_inst;
    float*                xf;            /* the instances' transforms, read by every draw */
    double                ctr[3], rad;   /* the bounding sphere of every instance */
    uint32_t*             colors;
    uint8_t*              sel;      /* the marks */
} scene_t;

/* "spacefill[:radius scale]", "licorice[:radius]", "ball[:ball scale]",
   "bonds[:radius]", "cartoon[:scale]", "ribbons[:half width]" */
static inline bool scene_parse_rep(const char* s, scene_rep_t* out) {
    char name[32] = {0};
    const char* colon = strchr(s, ':');
    size_t len = colon ? (size_t)(colon - s) : strlen(s);
    if (len >= sizeof(name)) return false;
    memcpy(name, s, len);
    memset(out, 0, sizeof(*out));
    out->param = colon ? (float)atof(colon + 1) : -1.0f;
    if      (!strcmp(name, "spacefill")) out->type = MD_RENDER_MOL_REP_SPACEFILL;
    else if (!strcmp(name, "licorice"))  out->type = MD_RENDER_MOL_REP_LICORICE;
    else if (!strcmp(name, "ball"))      out->type = MD_RENDER_MOL_REP_BALL_AND_STICK;
    else if (!strcmp(name, "bonds"))     out->type = MD_RENDER_MOL_REP_BONDS;
    else if (!strcmp(name, "cartoon"))   out->type = MD_RENDER_MOL_REP_CARTOON;
    else if (!strcmp(name, "ribbons"))   out->type = MD_RENDER_MOL_REP_RIBBONS;
    else { fprintf(stderr, "unknown representation '%s'\n", name); return false; }
    return true;
}

static inline bool scene_parse_fade_mode(const char* s, uint32_t* mode) {
    if (!strcmp(s, "thin"))   { *mode = MD_RENDER_FADE_THIN; return true; }
    if (!strcmp(s, "shrink")) { *mode = MD_RENDER_FADE_SHRINK; return true; }
    if (!strcmp(s, "both"))   { *mode = MD_RENDER_FADE_THIN | MD_RENDER_FADE_SHRINK; return true; }
    return false;
}

static inline void scene_mat_ident(float m[16]) { memset(m, 0, 64); m[0] = m[5] = m[10] = m[15] = 1.0f; }

/* A camera looking at `ctr` from direction `dir` (unit, from the target to the
   eye) at `dist`, reverse Z, infinite far plane, 45 degrees vertically */
static inline md_render_camera_t camera_look(const double ctr[3], const double dir[3], double dist, float aspect, double near_z) {
    md_render_camera_t c;
    memset(&c, 0, sizeof(c));
    double eye[3] = { ctr[0] + dir[0] * dist, ctr[1] + dir[1] * dist, ctr[2] + dir[2] * dist };
    double f[3] = { -dir[0], -dir[1], -dir[2] };                    /* forward */
    double up[3] = { 0, 1, 0 };
    if (fabs(f[1]) > 0.99) { up[1] = 0; up[2] = -1; }
    double r[3] = { f[1] * up[2] - f[2] * up[1], f[2] * up[0] - f[0] * up[2], f[0] * up[1] - f[1] * up[0] };
    double rl = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    for (int k = 0; k < 3; ++k) r[k] /= rl;
    double u[3] = { r[1] * f[2] - r[2] * f[1], r[2] * f[0] - r[0] * f[2], r[0] * f[1] - r[1] * f[0] };
    /* rows of the view rotation: r, u, -f (column-major storage) */
    c.view[0] = (float)r[0]; c.view[4] = (float)r[1]; c.view[8]  = (float)r[2];
    c.view[1] = (float)u[0]; c.view[5] = (float)u[1]; c.view[9]  = (float)u[2];
    c.view[2] = (float)-f[0]; c.view[6] = (float)-f[1]; c.view[10] = (float)-f[2];
    c.view[12] = (float)-(r[0] * eye[0] + r[1] * eye[1] + r[2] * eye[2]);
    c.view[13] = (float)-(u[0] * eye[0] + u[1] * eye[1] + u[2] * eye[2]);
    c.view[14] = (float)(f[0] * eye[0] + f[1] * eye[1] + f[2] * eye[2]);
    c.view[15] = 1.0f;
    const float tan_half = 0.4142f;                                 /* 45 degrees */
    c.proj[0] = 1.0f / (tan_half * aspect);
    c.proj[5] = 1.0f / tan_half;
    c.proj[11] = -1.0f;
    c.proj[14] = (float)near_z;
    return c;
}

/* The distance at which the bounding sphere just fits, and a near plane for it */
static inline double scene_fit_distance(const scene_t* sc, float aspect) {
    return sc->rad / sin(atan(0.4142 * (aspect < 1 ? aspect : 1.0)));
}
static inline double scene_near(const scene_t* sc, double dist) {
    return dist > sc->rad ? 0.05 * (dist - sc->rad) + 0.01 : 0.05;
}

/* The system, on the host: loaded, its coordinates replaced by `coords`',
   inferred */
static inline bool scene_load(scene_t* sc, const scene_desc_t* d) {
    memset(sc, 0, sizeof(*sc));
    sc->alloc = md_vm_arena_create(MD_GIGABYTES(16));
    sc->sys   = (md_system_t){ .alloc = sc->alloc };
    sc->state = (md_system_state_t){ .alloc = sc->alloc };
    md_str_t path = md_str_from_cstr(d->structure);
    const md_io_format_info_t* fmt = md_io_format_info(md_io_format_from_file(path));
    if (!fmt || !fmt->system_init || !fmt->system_init(&sc->sys, &sc->state, path, NULL)) {
        fprintf(stderr, "could not load '%s' as a system\n", d->structure);
        return false;
    }
    if (d->coords) {
        md_str_t cpath = md_str_from_cstr(d->coords);
        const md_io_format_info_t* cf = md_io_format_info(md_io_format_from_file(cpath));
        md_system_state_t cs = { .alloc = sc->alloc };
        if (!cf || !cf->state_init || !cf->state_init(&cs, cpath, sc->sys.atom.count) || cs.num_atoms != sc->sys.atom.count) {
            fprintf(stderr, "could not read coordinates for %zu atoms from '%s'\n", (size_t)sc->sys.atom.count, d->coords);
            return false;
        }
        sc->state = cs;
    }
    if (!md_system_infer(&sc->sys, &sc->state, NULL)) { fprintf(stderr, "md_system_infer failed\n"); return false; }
    sc->n = (uint32_t)sc->sys.atom.count;
    printf("%s: %u atoms, %zu bonds\n", d->structure, sc->n, (size_t)(sc->sys.bond.count + sc->sys.bond.interaction_count));
    return true;
}

/* Everything on the GPU, recorded in `s`: the topology, the atoms' colours
   and selection, the per-frame block (one regrouping update), the cartoon's
   profiles; then the batches of the representations */
static inline bool scene_upload(scene_t* sc, const scene_desc_t* d, md_gpu_stream_t s, md_render_mol_ctx_t kernels) {
    const uint32_t n = sc->n;
    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    if (!md_render_mol_plan(&plan, &sc->sys, NULL, temp)) { fprintf(stderr, "md_render_mol_plan failed\n"); md_temp_end(temp); return false; }
    if (!md_render_mol_write(s, &sc->topo, &plan, md_gpu_malloc(s, MD_GPU_MEM_DEVICE, plan.size).gpu, sc->content)) {
        fprintf(stderr, "md_render_mol_write failed\n");
        md_temp_end(temp);
        return false;
    }
    md_temp_end(temp);

    sc->colors = malloc((size_t)(n ? n : 1) * 4);
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t c = 0;
        if (sc->sys.atom.type_idx && sc->sys.atom.type.color && (size_t)sc->sys.atom.type_idx[i] < sc->sys.atom.type.count) c = sc->sys.atom.type.color[sc->sys.atom.type_idx[i]];
        if (!c) c = md_atomic_number_cpk_color(md_atom_atomic_number(&sc->sys.atom, i));
        sc->colors[i] = c | 0xFF000000u;
    }
    sc->sel = calloc(n ? n : 1, 1);
    uint32_t n_sel = 0;
    for (int64_t i = d->sel_beg < 0 ? 0 : d->sel_beg; i < d->sel_end && i < (int64_t)n; ++i) { sc->sel[i] = MD_RENDER_MARK_STATE0; n_sel++; }

    md_gpu_addr_t P = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, (size_t)n * 12).gpu;
    md_gpu_addr_t C = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, (size_t)n * 4).gpu;
    md_gpu_addr_t S = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, n ? n : 1).gpu;
    sc->dyn = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sc->topo.dynamic_size).gpu;
    if (!P || !C || !S || !sc->dyn) { fprintf(stderr, "out of GPU memory\n"); return false; }
    md_gpu_upload(s, P, sc->state.xyz, (size_t)n * 12);
    md_gpu_upload(s, C, sc->colors, (size_t)n * 4);
    md_gpu_upload(s, S, sc->sel, n ? n : 1);
    memset(&sc->atoms, 0, sizeof(sc->atoms));
    sc->atoms.position = P;
    md_render_mol_frame_cell(sc->atoms.cell, md_unitcell_flags(&sc->state.unitcell) ? &sc->state.unitcell : NULL);   /* the lipids whole across it */
    if (n_sel) sc->atoms.marks = S;
    if (!md_render_mol_update(kernels, s, &sc->topo, sc->dyn, &sc->atoms, MD_RENDER_MOL_UPDATE_REGROUP, sc->content)) { fprintf(stderr, "md_render_mol_update failed\n"); return false; }

    /* Instances along the cell (or the bounds) */
    const uint32_t replicate = d->replicate ? d->replicate : 1;
    double lo[3] = { 1e30, 1e30, 1e30 }, hi[3] = { -1e30, -1e30, -1e30 };
    for (uint32_t i = 0; i < n; ++i) {
        const float* p = (const float*)&sc->state.xyz[i];
        for (int k = 0; k < 3; ++k) { if (p[k] < lo[k]) lo[k] = p[k]; if (p[k] > hi[k]) hi[k] = p[k]; }
    }
    double A[3][3] = { { hi[0] - lo[0], 0, 0 }, { 0, hi[1] - lo[1], 0 }, { 0, 0, hi[2] - lo[2] } };
    if (md_unitcell_flags(&sc->state.unitcell) & (MD_UNITCELL_ORTHO | MD_UNITCELL_TRICLINIC)) md_unitcell_A_extract_double(A, &sc->state.unitcell);
    sc->n_inst = replicate * replicate * replicate;
    sc->xf = malloc((size_t)sc->n_inst * 64);
    for (uint32_t ix = 0, k = 0; ix < replicate; ++ix) for (uint32_t iy = 0; iy < replicate; ++iy) for (uint32_t iz = 0; iz < replicate; ++iz, ++k) {
        float* m = sc->xf + 16 * k;
        scene_mat_ident(m);
        for (int dd = 0; dd < 3; ++dd) m[12 + dd] = (float)(ix * A[0][dd] + iy * A[1][dd] + iz * A[2][dd]);
    }
    sc->rad = 0;
    for (int dd = 0; dd < 3; ++dd) {
        double ext = (replicate - 1) * (A[0][dd] + A[1][dd] + A[2][dd]);
        sc->ctr[dd] = 0.5 * (lo[dd] + hi[dd] + ext);
        double h = 0.5 * (hi[dd] - lo[dd] + fabs(ext));
        sc->rad += h * h;
    }
    sc->rad = sqrt(sc->rad);
    printf("%u instance(s), %llu atoms drawn; bounding radius %.1f A\n", sc->n_inst, (unsigned long long)n * sc->n_inst, sc->rad);

    /* The cartoon: the secondary structure (the state's, else computed), the
       profile per scale asked for, and its own colours: helix red, sheet
       yellow, coil grey (by the CA atoms); a nucleotide by its base, A red,
       C yellow, G green, T and U blue (all its atoms) */
    md_gpu_addr_t SSB = 0, CC = 0;
    if (sc->topo.backbone_count || sc->topo.nucleic_count) {
        const md_secondary_structure_t* ss = md_system_state_secondary_structure(&sc->state, &sc->sys);
        if (!ss && md_system_state_backbone_compute(&sc->state, &sc->sys)) ss = md_system_state_secondary_structure(&sc->state, &sc->sys);
        const size_t nseg = sc->sys.protein_backbone.segment.count;
        if (ss) {
            SSB = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, nseg * 4).gpu;
            md_gpu_upload(s, SSB, ss, nseg * 4);
        }
        uint32_t* cc = malloc((size_t)(n ? n : 1) * 4);
        for (uint32_t i = 0; i < n; ++i) cc[i] = sc->colors[i];
        for (size_t g = 0; g < nseg; ++g) {
            int32_t ca = sc->sys.protein_backbone.segment.atoms[g].ca;
            if (ca < 0 || (uint32_t)ca >= n) continue;
            int t = ss ? (int)ss[g] : 0;
            cc[ca] = (t >= MD_SECONDARY_STRUCTURE_HELIX_310 && t <= MD_SECONDARY_STRUCTURE_HELIX_PI) ? 0xFF4040E0u
                   : (t == MD_SECONDARY_STRUCTURE_BETA_SHEET || t == MD_SECONDARY_STRUCTURE_BETA_BRIDGE) ? 0xFF30C8F0u : 0xFFC8C8C8u;
        }
        const md_nucleic_backbone_data_t* NB = &sc->sys.nucleic_backbone;
        for (size_t g = 0; g < NB->segment.count; ++g) {
            const int32_t res = NB->segment.res_idx ? NB->segment.res_idx[g] : -1;
            if (res < 0) continue;
            md_str_t rn = md_system_residue_name(&sc->sys, (size_t)res);
            char base = 0;
            for (size_t k = 0; k < rn.len; ++k) if (strchr("ACGTU", rn.ptr[k])) base = rn.ptr[k];
            const uint32_t col = base == 'A' ? 0xFF4040E0u : base == 'C' ? 0xFF30C8F0u : base == 'G' ? 0xFF50C850u : base ? 0xFFE08040u : 0xFFC8C8C8u;
            md_urange_t ar = md_system_residue_atom_range(&sc->sys, (size_t)res);
            for (uint32_t i = (uint32_t)ar.beg; i < (uint32_t)ar.end && i < n; ++i) cc[i] = col;
        }
        CC = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, (size_t)n * 4).gpu;
        md_gpu_upload(s, CC, cc, (size_t)n * 4);
        free(cc);
        printf("cartoon: %u residues in %u pieces, secondary structure %s\n", sc->topo.backbone_count, sc->topo.backbone_chain_count, ss ? "given or computed" : "none (all coil)");
    }

    scene_rep_t def = { .type = MD_RENDER_MOL_REP_SPACEFILL, .param = -1.0f };
    const scene_rep_t* reps = d->rep_count ? d->reps : &def;
    const uint32_t n_reps = d->rep_count ? d->rep_count : 1;
    sc->batch_count = 0;
    for (uint32_t k = 0; k < n_reps; ++k) {
        md_render_mol_rep_t rep;
        memset(&rep, 0, sizeof(rep));
        rep.key  = k + 1;
        rep.type = reps[k].type;
        float p = reps[k].param;
        switch (rep.type) {
        case MD_RENDER_MOL_REP_SPACEFILL:
            rep.params.spacefill.radius_scale = p >= 0 ? p : 1.0f;
            memcpy(rep.params.spacefill.residues, reps[k].residues, sizeof(reps[k].residues));
            break;
        case MD_RENDER_MOL_REP_LICORICE:       rep.params.licorice.radius = p >= 0 ? p : 0.2f; break;
        case MD_RENDER_MOL_REP_BALL_AND_STICK: rep.params.ball_and_stick.ball_scale = p >= 0 ? p : 0.3f; rep.params.ball_and_stick.stick_radius = 0.12f; break;
        case MD_RENDER_MOL_REP_BONDS:          rep.params.bonds.radius = p >= 0 ? p : 0.15f; break;
        case MD_RENDER_MOL_REP_CARTOON: {
            const float sc_ = p >= 0 ? p : 1.0f, scale[3] = { sc_, sc_, sc_ };
            const size_t psize = md_render_mol_cartoon_size(&sc->topo);
            md_gpu_addr_t prof = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, psize ? psize : 256).gpu;
            if (sc->topo.backbone_count && !md_render_mol_cartoon_profile(kernels, s, &sc->topo, prof, SSB, scale, sc->content)) { fprintf(stderr, "md_render_mol_cartoon_profile failed\n"); return false; }
            rep.params.cartoon.profile = prof;
            break;
        }
        case MD_RENDER_MOL_REP_RIBBONS: rep.params.ribbons.half_width = p >= 0 ? p : 1.0f; break;
        default: break;
        }
        rep.color.mode  = MD_RENDER_COLOR_POINTS;
        rep.atom_color  = rep.type == MD_RENDER_MOL_REP_CARTOON && CC ? CC : C;
        rep.fade        = reps[k].fade;
        rep.chain_gaps  = d->chain_gaps;
        rep.bond_orders = d->bond_orders;
        if (sc->n_inst > 1) { rep.transform = sc->xf; rep.transform_count = sc->n_inst; }
        sc->batch_count += md_render_mol_emit(&sc->topo, sc->dyn, &sc->atoms, &rep, sc->batches + sc->batch_count, 64 - sc->batch_count);
    }
    return true;
}

static inline void scene_free(scene_t* sc) {
    free(sc->xf); free(sc->colors); free(sc->sel);
    md_system_free(&sc->sys);
    md_vm_arena_destroy(sc->alloc);
    memset(sc, 0, sizeof(*sc));
}

#endif

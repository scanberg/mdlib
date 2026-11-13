/*
md_render_bench: what clusters to cull by, and what md_render costs.

  md_render_bench clusters <structure> [<trajectory>] [--frames N] [--size N]
      CPU only. Loads a system (and a run of it), groups its atoms into
      clusters by several strategies, and reports how tightly each bounds its
      atoms, frame by frame. The measure is what culling pays for: the summed
      surface of the cluster bounding spheres per atom (a bounds test passes
      in proportion to it), and the largest radius (one huge cluster ruins
      occlusion around it).

  md_render_bench gpu [--atoms N[,N...]] [--res WxH] [--reps N] [--small P[,P...]] [--host-grid | --implicit] [--look NAME]
                     [--resident F | --resident-random F] [--mol-order cells [--jitter A]]
      Renders generated scenes (water-like boxes of three-atom molecules, in
      spacefill and licorice, near and far) through md_render, with each of
      the culling and raster options, and times full rasters, resolves and
      masks on the GPU stream, with the cull statistics beside them. By
      default the scene goes through md_render_mol (molecules in random
      order, clustered by md_render_mol_update, which is timed too, and
      the frame pipeline: a frame quantized on the host into GPU-visible
      memory, its upload into a slot of four with and without the regroup
      (and from host memory, for comparison), and the copy from 1, 2 or 4
      resident frames in a periodic box);
      --host-grid sorts and clusters by grid cell on the host instead, and
      --implicit sorts on the host and leaves md_render its runs of 64.
      --look keeps the looks whose name contains NAME ("licorice", "far").
      --resident F plans md_render_mol over that fraction of the molecules,
      one block in file order (as a structure is): what a view costs when the
      rest is not shown. --resident-random F spreads them at random instead,
      the worst case for the host's quantize (a set no representation makes).
      --mol-order cells generates the molecules in md_render_mol's grid-cell
      order instead of at random: the regroup's writes land in order.
      --jitter A then moves each molecule by up to A along each axis: the
      order of the last grouping, regrouped after the molecules moved.
      After each look's first row, one camera raster is profiled pass by
      pass (each pass synced, so each number carries a sync), with the cull
      by phase. MD_BENCH_SEQ=1 adds a sequence of rasters (a camera c0
      twice, c1 = c0 moved by MD_BENCH_SHIFT box widths (0.001), c1 again,
      back to c0) with the cull by phase for each: what a camera move costs
      phase 2.

Every run prints the mdlib revision it was built from, first.
*/

#include <md_system.h>
#include <md_system_state.h>
#include <md_hierarchy.h>
#include <md_infer.h>
#include <md_gro.h>
#include <md_pdb.h>
#include <md_xtc.h>
#include <md_trr.h>
#include <md_run.h>
#include <md_io.h>
#include <md_attributes.h>
#include <md_unitcell.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_str.h>

#if MD_ENABLE_GPU
#include <md_render.h>
#include <md_render_mol.h>
#include <core/md_gpu.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bench_rev.h"                    // MD_BENCH_SOURCE_REV, generated at build time by bench_version.cmake

#define ATOM_RADIUS 1.5                   // added to every cluster radius: a spacefill atom

/* =========================================================================
   Clusters
   ========================================================================= */

typedef struct clusters_t {
    uint32_t  count;
    uint32_t* offset;                     // [count + 1]
    uint32_t* atom;                       // [atoms]
} clusters_t;

typedef struct sys_info_t {
    const md_system_t* sys;
    size_t    n;
    int32_t*  mol_of_atom;                // molecule index, -1 none
    uint32_t* mol_size;                   // atoms per molecule
    uint32_t  small_max;                  // molecules of at most this many atoms diffuse
} sys_info_t;

static void clusters_alloc(clusters_t* c, size_t n) {
    c->count  = 0;
    c->offset = calloc(n + 1, sizeof(uint32_t));
    c->atom   = calloc(n, sizeof(uint32_t));
}

static void clusters_free(clusters_t* c) { free(c->offset); free(c->atom); memset(c, 0, sizeof(*c)); }

/* Greedy packing of units (runs of atoms) into clusters of at most S atoms; a
   unit larger than S is split into chunks of S. */
typedef struct packer_t {
    clusters_t* c;
    uint32_t    S;
    uint32_t    fill;      // atoms in the open cluster
    uint32_t    pos;       // atoms written
} packer_t;

static void pack_close(packer_t* p) {
    if (p->fill) { p->c->count++; p->c->offset[p->c->count] = p->pos; p->fill = 0; }
}

static void pack_unit(packer_t* p, const uint32_t* atoms, uint32_t n) {
    if (n > p->S) {
        pack_close(p);
        for (uint32_t i = 0; i < n; i += p->S) {
            uint32_t m = n - i < p->S ? n - i : p->S;
            memcpy(p->c->atom + p->pos, atoms + i, m * 4); p->pos += m; p->fill = m;
            pack_close(p);
        }
        return;
    }
    if (p->fill + n > p->S) pack_close(p);
    memcpy(p->c->atom + p->pos, atoms, n * 4);
    p->pos += n; p->fill += n;
}

/* A: chunks of S consecutive atoms */
static void build_chunks(clusters_t* c, const sys_info_t* si, uint32_t S) {
    clusters_alloc(c, si->n);
    for (size_t i = 0; i < si->n; ++i) c->atom[i] = (uint32_t)i;
    for (size_t i = 0; i < si->n; i += S) { c->count++; c->offset[c->count] = (uint32_t)(i + S < si->n ? i + S : si->n); }
}

/* B: whole residues, in file order, packed up to S, never across molecules */
static void build_residues(clusters_t* c, const sys_info_t* si, uint32_t S) {
    clusters_alloc(c, si->n);
    packer_t p = { c, S, 0, 0 };
    const md_residue_data_t* r = &si->sys->residue;
    uint32_t* tmp = malloc(si->n * 4);
    int32_t prev_mol = -2;
    if (r->count == 0) { build_chunks(c, si, S); free(tmp); return; }
    for (size_t k = 0; k < r->count; ++k) {
        uint32_t b = r->atom_offset[k], e = r->atom_offset[k + 1];
        int32_t mol = si->mol_of_atom[b];
        if (mol != prev_mol) pack_close(&p);
        prev_mol = mol;
        for (uint32_t i = b; i < e; ++i) tmp[i - b] = i;
        pack_unit(&p, tmp, e - b);
    }
    pack_close(&p);
    free(tmp);
}

/* Periodic helpers: the minimum image of d in a (possibly triclinic) cell */
typedef struct cell_t { bool on; double x, xy, xz, y, yz, z; } cell_t;

static cell_t cell_of(const md_unitcell_t* u) {
    cell_t c = {0};
    if (u && u->x > 0 && u->y > 0 && u->z > 0) {
        c.on = true; c.x = u->x; c.xy = u->xy; c.xz = u->xz; c.y = u->y; c.yz = u->yz; c.z = u->z;
    }
    return c;
}

static void min_image(double d[3], const cell_t* c) {
    if (!c->on) return;
    double fz = round(d[2] / c->z);
    d[0] -= fz * c->xz; d[1] -= fz * c->yz; d[2] -= fz * c->z;
    double fy = round(d[1] / c->y);
    d[0] -= fy * c->xy; d[1] -= fy * c->y;
    double fx = round(d[0] / c->x);
    d[0] -= fx * c->x;
}

static void unit_centroid(double out[3], const md_float3_t* xyz, const uint32_t* atoms, uint32_t n, const cell_t* cell) {
    double r[3] = { xyz[atoms[0]].x, xyz[atoms[0]].y, xyz[atoms[0]].z }, s[3] = {0};
    for (uint32_t i = 0; i < n; ++i) {
        double d[3] = { xyz[atoms[i]].x - r[0], xyz[atoms[i]].y - r[1], xyz[atoms[i]].z - r[2] };
        min_image(d, cell);
        s[0] += d[0]; s[1] += d[1]; s[2] += d[2];
    }
    out[0] = r[0] + s[0] / n; out[1] = r[1] + s[1] / n; out[2] = r[2] + s[2] / n;
    if (cell->on) {     /* into the box, so that Morton codes are of one image */
        double d[3] = { out[0] - 0.5 * cell->x, out[1] - 0.5 * cell->y, out[2] - 0.5 * cell->z };
        min_image(d, cell);
        out[0] = d[0] + 0.5 * cell->x; out[1] = d[1] + 0.5 * cell->y; out[2] = d[2] + 0.5 * cell->z;
    }
}

static uint32_t spread10(uint32_t x) {
    x &= 0x3FF;
    x = (x | (x << 16)) & 0x030000FF;
    x = (x | (x << 8))  & 0x0300F00F;
    x = (x | (x << 4))  & 0x030C30C3;
    x = (x | (x << 2))  & 0x09249249;
    return x;
}

typedef struct unit_t { uint32_t code; uint32_t begin, count; } unit_t;   // atoms in units_atoms[begin..)

static int unit_cmp(const void* a, const void* b) {
    uint32_t x = ((const unit_t*)a)->code, y = ((const unit_t*)b)->code;
    return x < y ? -1 : x > y;
}

/* Morton order of unit centroids, packed up to S. Units: the small molecules
   (only_small), or every residue of a large molecule and every small molecule. */
static void pack_morton(packer_t* p, const sys_info_t* si, const md_float3_t* xyz, const cell_t* cell, bool only_small) {
    const md_system_t* sys = si->sys;
    const md_molecule_data_t* m = &sys->molecule;
    uint32_t* ua = malloc(si->n * 4);
    unit_t* units = malloc((m->count + sys->residue.count + 1) * sizeof(unit_t));
    uint32_t nu = 0, na = 0;
    for (size_t k = 0; k < m->count; ++k) {
        uint32_t b = m->offset[k], e = m->offset[k + 1];
        if (e - b <= si->small_max) {
            units[nu].begin = na; units[nu].count = e - b;
            for (uint32_t s = b; s < e; ++s) ua[na++] = (uint32_t)m->atom_idx[s];
            nu++;
        }
    }
    if (!only_small) {
        const md_residue_data_t* r = &sys->residue;
        for (size_t k = 0; k < r->count; ++k) {
            uint32_t b = r->atom_offset[k], e = r->atom_offset[k + 1];
            int32_t mol = si->mol_of_atom[b];
            if (mol >= 0 && si->mol_size[mol] <= si->small_max) continue;
            units[nu].begin = na; units[nu].count = e - b;
            for (uint32_t i = b; i < e; ++i) ua[na++] = i;
            nu++;
        }
    }
    double lo[3] = { 1e30, 1e30, 1e30 }, hi[3] = { -1e30, -1e30, -1e30 };
    double* cen = malloc((nu + 1) * 3 * sizeof(double));
    for (uint32_t u = 0; u < nu; ++u) {
        unit_centroid(cen + 3 * u, xyz, ua + units[u].begin, units[u].count, cell);
        for (int k = 0; k < 3; ++k) { if (cen[3 * u + k] < lo[k]) lo[k] = cen[3 * u + k]; if (cen[3 * u + k] > hi[k]) hi[k] = cen[3 * u + k]; }
    }
    double ext = fmax(hi[0] - lo[0], fmax(hi[1] - lo[1], hi[2] - lo[2])) + 1e-6;
    for (uint32_t u = 0; u < nu; ++u) {
        uint32_t q[3];
        for (int k = 0; k < 3; ++k) q[k] = (uint32_t)((cen[3 * u + k] - lo[k]) / ext * 1023.0);
        units[u].code = spread10(q[0]) | (spread10(q[1]) << 1) | (spread10(q[2]) << 2);
    }
    qsort(units, nu, sizeof(unit_t), unit_cmp);
    for (uint32_t u = 0; u < nu; ++u) pack_unit(p, ua + units[u].begin, units[u].count);
    pack_close(p);
    free(cen); free(units); free(ua);
}

/* C: hybrid. Large molecules by residues in file order (as B, which stays
   valid over a trajectory: bonds keep it compact); small molecules by Morton
   order of their centroids, at the frame given (rebuilt per frame or not). */
static void build_hybrid(clusters_t* c, const sys_info_t* si, uint32_t S, const md_float3_t* xyz, const cell_t* cell) {
    clusters_alloc(c, si->n);
    packer_t p = { c, S, 0, 0 };
    const md_residue_data_t* r = &si->sys->residue;
    uint32_t* tmp = malloc(si->n * 4);
    int32_t prev_mol = -2;
    for (size_t k = 0; k < r->count; ++k) {
        uint32_t b = r->atom_offset[k], e = r->atom_offset[k + 1];
        int32_t mol = si->mol_of_atom[b];
        if (mol >= 0 && si->mol_size[mol] <= si->small_max) continue;
        if (mol != prev_mol) pack_close(&p);
        prev_mol = mol;
        for (uint32_t i = b; i < e; ++i) tmp[i - b] = i;
        pack_unit(&p, tmp, e - b);
    }
    pack_close(&p);
    free(tmp);
    pack_morton(&p, si, xyz, cell, true);
    /* atoms in no molecule and no residue would be lost: none in practice; checked by the caller */
}

/* E: hybrid with a uniform grid for the small molecules: cells sized for
   about 3/4 of S atoms at the system's density, one cluster per cell (packed
   up to S, so a full cell spills into a second). What a GPU counting sort
   builds in O(n). */
static void build_grid(clusters_t* c, const sys_info_t* si, uint32_t S, const md_float3_t* xyz, const cell_t* cell) {
    clusters_alloc(c, si->n);
    packer_t p = { c, S, 0, 0 };
    const md_system_t* sys = si->sys;
    /* large molecules as in hybrid */
    const md_residue_data_t* r = &sys->residue;
    uint32_t* tmp = malloc(si->n * 4);
    int32_t prev_mol = -2;
    for (size_t k = 0; k < r->count; ++k) {
        uint32_t b = r->atom_offset[k], e = r->atom_offset[k + 1];
        int32_t mol = si->mol_of_atom[b];
        if (mol >= 0 && si->mol_size[mol] <= si->small_max) continue;
        if (mol != prev_mol) pack_close(&p);
        prev_mol = mol;
        for (uint32_t i = b; i < e; ++i) tmp[i - b] = i;
        pack_unit(&p, tmp, e - b);
    }
    pack_close(&p);
    free(tmp);
    /* small molecules: centroids, bounds, cells */
    const md_molecule_data_t* m = &sys->molecule;
    uint32_t nu = 0, small_atoms = 0;
    for (size_t k = 0; k < m->count; ++k) if (m->offset[k + 1] - m->offset[k] <= si->small_max) { nu++; small_atoms += m->offset[k + 1] - m->offset[k]; }
    if (!nu) return;
    uint32_t* mol = malloc(nu * 4);
    double* cen = malloc(nu * 3 * sizeof(double));
    uint32_t* ua = malloc(si->n * 4);
    double lo[3] = { 1e30, 1e30, 1e30 }, hi[3] = { -1e30, -1e30, -1e30 };
    for (size_t k = 0, u = 0; k < m->count; ++k) {
        uint32_t b = m->offset[k], e = m->offset[k + 1];
        if (e - b > si->small_max) continue;
        for (uint32_t s2 = b; s2 < e; ++s2) ua[s2 - b] = (uint32_t)m->atom_idx[s2];
        unit_centroid(cen + 3 * u, xyz, ua, e - b, cell);
        for (int j = 0; j < 3; ++j) { if (cen[3 * u + j] < lo[j]) lo[j] = cen[3 * u + j]; if (cen[3 * u + j] > hi[j]) hi[j] = cen[3 * u + j]; }
        mol[u++] = (uint32_t)k;
    }
    double vol = (hi[0] - lo[0] + 1e-3) * (hi[1] - lo[1] + 1e-3) * (hi[2] - lo[2] + 1e-3);
    double edge = cbrt(0.75 * S * vol / small_atoms);
    uint32_t dim[3];
    for (int j = 0; j < 3; ++j) { dim[j] = (uint32_t)((hi[j] - lo[j]) / edge) + 1; }
    size_t ncell = (size_t)dim[0] * dim[1] * dim[2];
    uint32_t* cnt = calloc(ncell + 1, 4);
    uint32_t* cell_of_u = malloc(nu * 4);
    for (uint32_t u = 0; u < nu; ++u) {
        uint32_t q[3];
        for (int j = 0; j < 3; ++j) { q[j] = (uint32_t)((cen[3 * u + j] - lo[j]) / edge); if (q[j] >= dim[j]) q[j] = dim[j] - 1; }
        cell_of_u[u] = (q[2] * dim[1] + q[1]) * dim[0] + q[0];
        cnt[cell_of_u[u] + 1]++;
    }
    for (size_t i = 0; i < ncell; ++i) cnt[i + 1] += cnt[i];
    uint32_t* order = malloc(nu * 4);
    uint32_t* fill = calloc(ncell, 4);
    for (uint32_t u = 0; u < nu; ++u) order[cnt[cell_of_u[u]] + fill[cell_of_u[u]]++] = u;
    for (size_t ci = 0; ci < ncell; ++ci) {
        for (uint32_t j = cnt[ci]; j < cnt[ci + 1]; ++j) {
            uint32_t k = mol[order[j]], b = m->offset[k], e = m->offset[k + 1];
            for (uint32_t s2 = b; s2 < e; ++s2) ua[s2 - b] = (uint32_t)m->atom_idx[s2];
            pack_unit(&p, ua, e - b);
        }
        pack_close(&p);
    }
    free(mol); free(cen); free(ua); free(cnt); free(cell_of_u); free(order); free(fill);
}

/* D: everything spatial (residues of large molecules and small molecules, Morton order) */
static void build_morton(clusters_t* c, const sys_info_t* si, uint32_t S, const md_float3_t* xyz, const cell_t* cell) {
    clusters_alloc(c, si->n);
    packer_t p = { c, S, 0, 0 };
    pack_morton(&p, si, xyz, cell, false);
}

typedef struct quality_t {
    double area_per_atom;    // sum of 4 pi r^2 over clusters / atoms
    double mean_r, p95_r, max_r;
    uint32_t clusters;
    uint32_t atoms;
} quality_t;

static int dcmp(const void* a, const void* b) { double x = *(const double*)a, y = *(const double*)b; return x < y ? -1 : x > y; }

static quality_t evaluate(const clusters_t* c, const md_float3_t* xyz, const cell_t* cell) {
    quality_t q = {0};
    double* rad = malloc((c->count + 1) * sizeof(double));
    double area = 0, sum = 0;
    for (uint32_t k = 0; k < c->count; ++k) {
        uint32_t b = c->offset[k], e = c->offset[k + 1];
        double cen[3];
        /* the centroid of the cluster's atoms, each at its image nearest the first */
        double r0[3] = { xyz[c->atom[b]].x, xyz[c->atom[b]].y, xyz[c->atom[b]].z }, s[3] = {0};
        for (uint32_t i = b; i < e; ++i) {
            double d[3] = { xyz[c->atom[i]].x - r0[0], xyz[c->atom[i]].y - r0[1], xyz[c->atom[i]].z - r0[2] };
            min_image(d, cell);
            s[0] += d[0]; s[1] += d[1]; s[2] += d[2];
        }
        for (int j = 0; j < 3; ++j) cen[j] = r0[j] + s[j] / (e - b);
        double r = 0;
        for (uint32_t i = b; i < e; ++i) {
            double d[3] = { xyz[c->atom[i]].x - cen[0], xyz[c->atom[i]].y - cen[1], xyz[c->atom[i]].z - cen[2] };
            min_image(d, cell);
            double l = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            if (l > r) r = l;
        }
        r += ATOM_RADIUS;
        rad[k] = r;
        area += 4.0 * MD_PI * r * r;
        sum += r;
        q.atoms += e - b;
    }
    qsort(rad, c->count, sizeof(double), dcmp);
    q.clusters = c->count;
    q.area_per_atom = q.atoms ? area / q.atoms : 0;
    q.mean_r = c->count ? sum / c->count : 0;
    q.p95_r  = c->count ? rad[(size_t)(0.95 * (c->count - 1))] : 0;
    q.max_r  = c->count ? rad[c->count - 1] : 0;
    free(rad);
    return q;
}

static bool has_ext(const char* path, const char* ext) {
    size_t n = strlen(path), m = strlen(ext);
    return n >= m && strcmp(path + n - m, ext) == 0;
}

static int run_clusters(int argc, char** argv) {
    const char* structure = NULL;
    const char* traj = NULL;
    uint32_t max_frames = 40;
    uint32_t sizes[3] = { 32, 64, 128 };
    uint32_t nsizes = 3;
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i], "--frames") && i + 1 < argc) max_frames = (uint32_t)atoi(argv[++i]);
        else if (!strcmp(argv[i], "--size") && i + 1 < argc) { sizes[0] = (uint32_t)atoi(argv[++i]); nsizes = 1; }
        else if (!structure) structure = argv[i];
        else traj = argv[i];
    }
    if (!structure) { fprintf(stderr, "usage: md_render_bench clusters <structure> [<trajectory>] [--frames N] [--size N]\n"); return 1; }

    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(4));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    bool ok = false;
    md_str_t path = md_str_from_cstr(structure);
    if (has_ext(structure, ".gro"))      ok = md_gro_system_init_from_file(&sys, &state, path, NULL);
    else if (has_ext(structure, ".pdb")) ok = md_pdb_system_init_from_file(&sys, &state, path, NULL);
    if (!ok) { fprintf(stderr, "could not load %s (gro or pdb)\n", structure); return 1; }
    md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL);

    const md_str_t run = MD_STR_LIT("bench");
    size_t num_frames = 1;
    md_system_extract_t* ex = NULL;
    if (traj) {
        md_str_t tp = md_str_from_cstr(traj);
        bool pub = has_ext(traj, ".xtc") ? md_xtc_system_publish_run(&sys, tp, run, MD_RUN_FLAG_DISABLE_CACHE_WRITE)
                 : has_ext(traj, ".trr") ? md_trr_system_publish_run(&sys, tp, run, MD_RUN_FLAG_DISABLE_CACHE_WRITE)
                 : has_ext(traj, ".pdb") ? md_pdb_system_publish_run(&sys, tp, run, MD_RUN_FLAG_DISABLE_CACHE_WRITE) : false;
        if (!pub) { fprintf(stderr, "could not open the run %s\n", traj); return 1; }
        char buf[256];
        const md_attribute_t* t = md_attributes_find(&sys.attributes, md_run_path(buf, sizeof(buf), run, MD_STR_LIT("time")));
        num_frames = t ? t->format.shape[0] : 0;
        const md_str_t paths[] = { MD_STR_LIT("atom/position"), MD_STR_LIT("unitcell") };
        ex = md_system_extract_begin(&sys, run, paths, 2, md_get_heap_allocator());
        if (!ex || !num_frames) { fprintf(stderr, "no frames in %s\n", traj); return 1; }
    }

    sys_info_t si = { .sys = &sys, .n = sys.atom.count, .small_max = 16 };
    si.mol_of_atom = malloc(si.n * sizeof(int32_t));
    si.mol_size    = calloc(sys.molecule.count + 1, sizeof(uint32_t));
    for (size_t i = 0; i < si.n; ++i) si.mol_of_atom[i] = -1;
    size_t small_atoms = 0;
    for (size_t k = 0; k < sys.molecule.count; ++k) {
        uint32_t b = sys.molecule.offset[k], e = sys.molecule.offset[k + 1];
        si.mol_size[k] = e - b;
        for (uint32_t s = b; s < e; ++s) si.mol_of_atom[sys.molecule.atom_idx[s]] = (int32_t)k;
        if (e - b <= si.small_max) small_atoms += e - b;
    }

    printf("mdlib %s\n", MD_BENCH_SOURCE_REV);
    printf("%s%s%s: %zu atoms, %zu residues, %zu molecules (%zu atoms in molecules of <= %u atoms), %zu frames%s\n",
           structure, traj ? " + " : "", traj ? traj : "", si.n, sys.residue.count, sys.molecule.count, small_atoms, si.small_max,
           num_frames, num_frames > max_frames ? " (sampled)" : "");
    printf("measure: area = sum over clusters of 4 pi r^2 per atom (A^2), r = bounding sphere + %.1f A; lower is better\n\n", ATOM_RADIUS);

    enum { A_CHUNKS, B_RESIDUES, C_HYBRID_STATIC, C_HYBRID, E_GRID, D_MORTON, NSTRAT };
    const char* names[NSTRAT] = { "chunks (static)", "residues (static)", "hybrid, built once", "hybrid, per frame", "grid, per frame", "morton, per frame" };

    uint32_t nf = num_frames < max_frames ? (uint32_t)num_frames : max_frames;
    for (uint32_t si_ = 0; si_ < nsizes; ++si_) {
        const uint32_t S = sizes[si_];
        printf("cluster size <= %u\n", S);
        printf("  %-20s %9s | %8s %8s %8s %8s | %8s %8s %8s | %9s\n", "strategy", "clusters", "area f0", "mean r", "p95 r", "max r",
               "area end", "p95 end", "max end", "worst area");
        clusters_t st[NSTRAT] = {0};
        quality_t first[NSTRAT], last[NSTRAT];
        double worst[NSTRAT] = {0};
        for (uint32_t fi = 0; fi < nf; ++fi) {
            int64_t frame = nf > 1 ? (int64_t)((double)fi * (num_frames - 1) / (nf - 1)) : 0;
            if (ex && !md_system_extract_frame(ex, frame, &state)) { fprintf(stderr, "frame %lld failed\n", (long long)frame); return 1; }
            cell_t cell = cell_of(&state.unitcell);
            for (int s = 0; s < NSTRAT; ++s) {
                bool rebuild = fi == 0 || s == C_HYBRID || s == D_MORTON || s == E_GRID;
                if (rebuild) {
                    clusters_free(&st[s]);
                    switch (s) {
                    case A_CHUNKS:        build_chunks(&st[s], &si, S); break;
                    case B_RESIDUES:      build_residues(&st[s], &si, S); break;
                    case C_HYBRID_STATIC:
                    case C_HYBRID:        build_hybrid(&st[s], &si, S, state.xyz, &cell); break;
                    case D_MORTON:        build_morton(&st[s], &si, S, state.xyz, &cell); break;
                    case E_GRID:          build_grid(&st[s], &si, S, state.xyz, &cell); break;
                    }
                }
                quality_t q = evaluate(&st[s], state.xyz, &cell);
                if (q.atoms != si.n && fi == 0) printf("  (%s covers %u of %zu atoms)\n", names[s], q.atoms, si.n);
                if (fi == 0) first[s] = q;
                last[s] = q;
                if (q.area_per_atom > worst[s]) worst[s] = q.area_per_atom;
            }
        }
        for (int s = 0; s < NSTRAT; ++s) {
            printf("  %-20s %9u | %8.1f %8.2f %8.2f %8.2f | %8.1f %8.2f %8.2f | %9.1f\n", names[s], first[s].clusters,
                   first[s].area_per_atom, first[s].mean_r, first[s].p95_r, first[s].max_r,
                   last[s].area_per_atom, last[s].p95_r, last[s].max_r, worst[s]);
            clusters_free(&st[s]);
        }
        printf("\n");
    }
    if (ex) md_system_extract_end(ex);
    free(si.mol_of_atom); free(si.mol_size);
    md_system_free(&sys);
    return 0;
}

/* =========================================================================
   GPU timings
   ========================================================================= */

#if MD_ENABLE_GPU

typedef enum cluster_src_t {
    SRC_MOL,          /* md_render_mol: the molecules in random order, regrouped by md_render_mol_update */
    SRC_HOST_GRID,    /* sorted by grid cell on the host, a cluster per cell */
    SRC_IMPLICIT,     /* sorted by grid cell on the host, runs of 64 */
} cluster_src_t;

typedef struct scene_t {
    uint32_t      atoms, bonds;
    md_gpu_addr_t pos, rad, bond_idx, color, sel;
    uint32_t      clusters;          /* SRC_HOST_GRID: grid cells */
    md_gpu_addr_t atom_clusters, bond_clusters;
    float         box;            // edge, A
    /* SRC_MOL */
    md_render_mol_topo_t topo;
    md_gpu_addr_t        dyn;            /* with room for velocities */
    float*               host_pos;       /* the positions, for the frames */
    md_bitmask_t         resident;       /* --resident F < 1: the resident atoms */
    bool                 partial;
    double               plan_ms, write_ms;  /* what a new plan costs (a change of the resident set) */
    md_render_mol_plan_t plan;           /* a copy of the plan: written again, the first regroup is again the first */
} scene_t;

typedef struct mol_cell_t { uint32_t cell, m; } mol_cell_t;
static int mol_cell_cmp(const void* a, const void* b) {
    uint32_t x = ((const mol_cell_t*)a)->cell, y = ((const mol_cell_t*)b)->cell;
    return x < y ? -1 : x > y;
}

/* A box of three-atom molecules at the density of water (0.1 atoms / A^3):
   O anywhere, two H at 1 A from it. The molecules are ordered by the grid
   cell of their O, cells sized for about 64 atoms, and each cell is a
   cluster: what md_render_mol builds for solvent. */
/* --mol-order cells: the molecules sorted by the cell of their first atom,
   z, then y, then x (md_render_mol's grid numbering), at about its cell size */
typedef struct mol_key_t { uint64_t key; uint32_t m; } mol_key_t;
static int cmp_mol_key(const void* a, const void* b) {
    const mol_key_t* x = a; const mol_key_t* y = b;
    return x->key < y->key ? -1 : x->key > y->key ? 1 : (x->m < y->m ? -1 : x->m > y->m);
}

static bool scene_make(md_gpu_stream_t s, scene_t* sc, uint32_t molecules, cluster_src_t src, float resident, bool scattered, bool cell_order, float jitter) {
    memset(sc, 0, sizeof(*sc));
    sc->atoms = molecules * 3;
    sc->bonds = molecules * 2;
    sc->box   = (float)cbrt(sc->atoms / 0.1);
    float* pos = malloc((size_t)sc->atoms * 12);
    float* rad = malloc((size_t)sc->atoms * 4);
    uint32_t* idx = malloc((size_t)sc->bonds * 8);
    uint32_t* col = malloc((size_t)sc->atoms * 4);
    uint32_t seed = 7;
    #define RND() (((seed = seed * 1664525u + 1013904223u) >> 8) / 16777216.0f)
    for (uint32_t m = 0; m < molecules; ++m) {
        float o[3] = { RND() * sc->box, RND() * sc->box, RND() * sc->box };
        for (int a = 0; a < 3; ++a) {
            uint32_t i = 3 * m + a;
            float d[3] = { 0, 0, 0 };
            if (a) { d[0] = RND() - 0.5f; d[1] = RND() - 0.5f; d[2] = RND() - 0.5f;
                     float l = sqrtf(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) + 1e-6f; d[0] /= l; d[1] /= l; d[2] /= l; }
            for (int k = 0; k < 3; ++k) pos[3 * i + k] = o[k] + d[k];
            rad[i] = a ? 1.1f : 1.52f;
            col[i] = a ? 0xFFFFFFFFu : 0xFF0000FFu;
        }
        idx[4 * m + 0] = 3 * m; idx[4 * m + 1] = 3 * m + 1;
        idx[4 * m + 2] = 3 * m; idx[4 * m + 3] = 3 * m + 2;
    }
    #undef RND
    if (cell_order) {
        const float e = (float)cbrt(64.0 / 0.1);           /* CELL_TARGET atoms at 0.1 per A^3 */
        const uint64_t d = (uint64_t)(sc->box / e) + 1;
        mol_key_t* k = malloc((size_t)molecules * sizeof(mol_key_t));
        for (uint32_t m = 0; m < molecules; ++m) {
            const float* o = pos + 9 * (size_t)m;          /* atom 0 sits at the origin */
            uint64_t c[3];
            for (int j = 0; j < 3; ++j) { float q = floorf(o[j] / e); c[j] = q < 0 ? 0 : (uint64_t)q; }
            k[m].key = (c[2] * d + c[1]) * d + c[0];
            k[m].m = m;
        }
        qsort(k, molecules, sizeof(mol_key_t), cmp_mol_key);
        float* p2 = malloc((size_t)sc->atoms * 12);
        float* r2 = malloc((size_t)sc->atoms * 4);
        uint32_t* c2 = malloc((size_t)sc->atoms * 4);
        for (uint32_t m = 0; m < molecules; ++m) {
            const uint32_t f = k[m].m;
            memcpy(p2 + 9 * (size_t)m, pos + 9 * (size_t)f, 36);
            memcpy(r2 + 3 * (size_t)m, rad + 3 * (size_t)f, 12);
            memcpy(c2 + 3 * (size_t)m, col + 3 * (size_t)f, 12);
        }
        free(pos); free(rad); free(col); free(k);
        pos = p2; rad = r2; col = c2;
        if (jitter > 0.0f) {
            /* --jitter A: then each molecule moved by up to A along each axis,
               after the sort: the order of a grouping, molecules since moved
               (what a regroup sees at the next) */
            uint32_t js = 12345;
            #define JRND() (((js = js * 1664525u + 1013904223u) >> 8) / 16777216.0f)
            for (uint32_t m = 0; m < molecules; ++m) {
                const float d[3] = { (2.0f * JRND() - 1.0f) * jitter, (2.0f * JRND() - 1.0f) * jitter, (2.0f * JRND() - 1.0f) * jitter };
                for (int a = 0; a < 3; ++a) for (int j = 0; j < 3; ++j) pos[9 * (size_t)m + 3 * a + j] += d[j];
            }
            #undef JRND
        }
    }
    if (src == SRC_MOL) {
        /* A system of waters for md_render_mol: molecules and bonds, nothing else */
        md_system_t sys;
        memset(&sys, 0, sizeof(sys));
        sys.topology_version = 1;
        sys.atom.count = sc->atoms;
        uint32_t* moff = malloc(((size_t)molecules + 1) * 4);
        int32_t*  mats = malloc((size_t)sc->atoms * 4);
        md_atom_pair_t* pairs = malloc((size_t)sc->bonds * sizeof(md_atom_pair_t));
        for (uint32_t m = 0; m <= molecules; ++m) moff[m] = 3 * m;
        for (uint32_t i = 0; i < sc->atoms; ++i) mats[i] = (int32_t)i;
        for (uint32_t b = 0; b < sc->bonds; ++b) { pairs[b].idx[0] = (md_atom_idx_t)idx[2 * b]; pairs[b].idx[1] = (md_atom_idx_t)idx[2 * b + 1]; }
        sys.molecule.count = molecules; sys.molecule.offset = moff; sys.molecule.atom_idx = mats;
        /* two particle types, O and H: the radii come from their table */
        md_particle_type_idx_t* types = malloc((size_t)sc->atoms * sizeof(md_particle_type_idx_t));
        for (uint32_t i = 0; i < sc->atoms; ++i) types[i] = (i % 3) ? 2 : 1;
        float type_radius[3] = { 0.0f, 1.52f, 1.1f };
        sys.atom.type_idx = types;
        sys.atom.type.count = 3; sys.atom.type.radius = type_radius;
        sys.bond.count = sc->bonds; sys.bond.pairs = pairs;
        /* --resident: that fraction of the molecules, one block in file order
           (as a structure is); --resident-random: spread at random (the worst
           case, a set no representation makes) */
        md_bitset_t rset;
        if (resident < 1.0f) {
            sc->partial  = true;
            sc->resident = md_bitmask_create(0, sc->atoms, md_get_heap_allocator());
            if (scattered) {
                const uint32_t cut = (uint32_t)(resident * 16777216.0f);
                for (uint32_t m = 0; m < molecules; ++m) {
                    if (((m * 2654435761u) >> 8) >= cut) continue;
                    md_bitmask_set_range(&sc->resident, 3 * m, 3 * m + 3);
                }
            } else {
                const uint32_t mb = (uint32_t)(molecules * (1.0f - resident) * 0.5f);
                md_bitmask_set_range(&sc->resident, 3 * mb, 3 * (mb + (uint32_t)(molecules * resident)));
            }
            rset = md_bitmask_view(&sc->resident);
        }
        md_temp_scope_t temp = md_temp_begin();
        md_render_mol_plan_t plan;
        md_tick_t t0 = md_tick_now();
        bool ok = md_render_mol_plan(&plan, &sys, sc->partial ? &rset : NULL, temp);
        md_tick_t t1 = md_tick_now();
        md_gpu_addr_t base = ok ? md_gpu_malloc(s, MD_GPU_MEM_DEVICE, plan.size).gpu : 0;
        ok = ok && base && md_render_mol_write(s, &sc->topo, &plan, base, NULL);
        md_gpu_stream_sync(s);
        sc->plan_ms  = md_tick_to_milliseconds(t1 - t0);
        sc->write_ms = md_tick_to_milliseconds(md_tick_now() - t1);
        if (ok) {
            void* keep = malloc(plan.size);
            if (keep) { memcpy(keep, plan.data, plan.size); sc->plan = plan; sc->plan.data = keep; }
        }
        md_temp_end(temp);
        sc->dyn = ok ? md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sc->topo.dynamic_size_velocity).gpu : 0;
        free(moff); free(mats); free(pairs); free(types);
        if (!ok || !sc->dyn) { free(pos); free(rad); free(idx); free(col); return false; }
    } else {
    /* order by cell (Morton order of the cells), and the clusters */
        const float edge = cbrtf(64.0f / 0.1f);
        uint32_t dim = (uint32_t)(sc->box / edge) + 1;
        mol_cell_t* mc = malloc((size_t)molecules * sizeof(mol_cell_t));
        for (uint32_t m = 0; m < molecules; ++m) {
            uint32_t q[3];
            for (int k = 0; k < 3; ++k) { q[k] = (uint32_t)(pos[9 * m + k] / edge); if (q[k] >= dim) q[k] = dim - 1; }
            uint32_t code = 0;
            for (int bit = 0; bit < 10; ++bit) for (int k = 0; k < 3; ++k) code |= ((q[k] >> bit) & 1u) << (3 * bit + k);
            mc[m].cell = code; mc[m].m = m;
        }
        qsort(mc, molecules, sizeof(mol_cell_t), mol_cell_cmp);
        float* p2 = malloc((size_t)sc->atoms * 12);
        float* r2 = malloc((size_t)sc->atoms * 4);
        uint32_t* c2 = malloc((size_t)sc->atoms * 4);
        uint32_t* aoff = malloc(((size_t)molecules + 1) * 4);
        uint32_t* boff = malloc(((size_t)molecules + 1) * 4);
        uint32_t nc = 0;
        for (uint32_t j = 0; j < molecules; ++j) {
            uint32_t m = mc[j].m;
            if (j == 0 || mc[j].cell != mc[j - 1].cell) { aoff[nc] = 3 * j; boff[nc] = 2 * j; nc++; }
            memcpy(p2 + 9 * j, pos + 9 * m, 36);
            memcpy(r2 + 3 * j, rad + 3 * m, 12);
            memcpy(c2 + 3 * j, col + 3 * m, 12);
        }
        aoff[nc] = sc->atoms; boff[nc] = sc->bonds;
        memcpy(pos, p2, (size_t)sc->atoms * 12); memcpy(rad, r2, (size_t)sc->atoms * 4); memcpy(col, c2, (size_t)sc->atoms * 4);
        sc->clusters = nc;
        sc->atom_clusters = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, ((size_t)nc + 1) * 4).gpu;
        sc->bond_clusters = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, ((size_t)nc + 1) * 4).gpu;
        md_gpu_upload(s, sc->atom_clusters, aoff, ((size_t)nc + 1) * 4);
        md_gpu_upload(s, sc->bond_clusters, boff, ((size_t)nc + 1) * 4);
        free(mc); free(p2); free(r2); free(c2); free(aoff); free(boff);
    }
    size_t words = (sc->atoms + 31) / 32;
    uint32_t* sel = calloc(words, 4);
    for (uint32_t i = 0; i < sc->atoms; i += 97) sel[i / 32] |= 1u << (i % 32);     // about 1 % selected
    sc->pos      = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, (size_t)sc->atoms * 12).gpu;
    /* md_render_mol reads radii from its types and bonds from its topology:
       the scene's own arrays only for the plain batches */
    const bool plain = src != SRC_MOL;
    sc->rad      = plain ? md_gpu_malloc(s, MD_GPU_MEM_DEVICE, (size_t)sc->atoms * 4).gpu : 0;
    sc->bond_idx = plain ? md_gpu_malloc(s, MD_GPU_MEM_DEVICE, (size_t)sc->bonds * 8).gpu : 0;
    sc->color    = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, (size_t)sc->atoms * 4).gpu;
    sc->sel      = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, words * 4).gpu;
    bool ok = sc->pos && (!plain || (sc->rad && sc->bond_idx)) && sc->color && sc->sel &&
              md_gpu_upload(s, sc->pos, pos, (size_t)sc->atoms * 12) &&
              (!plain || (md_gpu_upload(s, sc->rad, rad, (size_t)sc->atoms * 4) && md_gpu_upload(s, sc->bond_idx, idx, (size_t)sc->bonds * 8))) &&
              md_gpu_upload(s, sc->color, col, (size_t)sc->atoms * 4) && md_gpu_upload(s, sc->sel, sel, words * 4);
    md_gpu_stream_sync(s);
    if (src == SRC_MOL) sc->host_pos = pos; else free(pos);
    free(rad); free(idx); free(col); free(sel);
    return ok;
}

static void scene_free(md_gpu_stream_t s, scene_t* sc) {
    if (sc->atom_clusters) md_gpu_free(s, sc->atom_clusters);
    if (sc->bond_clusters) md_gpu_free(s, sc->bond_clusters);
    if (sc->topo.base) md_gpu_free(s, sc->topo.base);
    if (sc->dyn) md_gpu_free(s, sc->dyn);
    md_gpu_free(s, sc->pos); md_gpu_free(s, sc->color); md_gpu_free(s, sc->sel);
    if (sc->rad) md_gpu_free(s, sc->rad);
    if (sc->bond_idx) md_gpu_free(s, sc->bond_idx);
    free(sc->host_pos);
    free((void*)sc->plan.data);
    if (sc->partial) md_bitmask_free(&sc->resident);
}

static md_render_camera_t camera_at(const scene_t* sc, float distance_in_boxes, float aspect) {
    md_render_camera_t c;
    memset(&c, 0, sizeof(c));
    const float tan_half = 0.4142f;    /* 45 degrees */
    /* look down -Z at the box centre from in front of its +Z face */
    float cx = 0.5f * sc->box, cy = 0.5f * sc->box, cz = sc->box + distance_in_boxes * sc->box;
    c.view[0] = c.view[5] = c.view[10] = c.view[15] = 1.0f;
    c.view[12] = -cx; c.view[13] = -cy; c.view[14] = -cz;
    c.proj[0] = 1.0f / (tan_half * aspect);
    c.proj[5] = 1.0f / tan_half;
    c.proj[11] = -1.0f;
    c.proj[14] = 0.1f;
    return c;
}

static int cmp_d(const void* a, const void* b) { double x = *(const double*)a, y = *(const double*)b; return x < y ? -1 : x > y; }

/* md_render_mol_update with a hook after each launch: in md_render_mol.c,
   not in its header (it is for this profile only) */
bool md_render_mol_update_hooked(md_render_mol_t mol, md_gpu_stream_t stream, const md_render_mol_topo_t* topo,
                                 md_gpu_addr_t dynamic, const md_render_mol_atoms_t* atoms,
                                 md_render_mol_update_flags_t flags, void (*hook)(void* user, const char* pass), void* user);

bool md_render_mol_frame_upload_hooked(md_render_mol_t mol, md_gpu_stream_t stream, const md_render_mol_topo_t* topo, md_gpu_addr_t dynamic,
                                       md_render_mol_frames_t* frames, uint32_t slot, const md_render_mol_frame_t* frame,
                                       md_render_mol_update_flags_t flags, void (*hook)(void* user, const char* pass), void* user);

/* md_render's per-pass hook and the phase-1 occlusion margin experiment: in
   md_render.c, not in its header */
void md_render_debug_pass_hook(void (*hook)(void* user, const char* pass), void* user);
void md_render_mol_debug_prefetch(uint32_t windows_ahead);
void md_render_mol_debug_write(uint32_t leave_out);
void md_render_mol_debug_frame_groups(uint32_t groups);
bool md_render_debug_kernel_info(md_render_t r, uint32_t i, const char** name, md_gpu_kernel_info_t* info);
bool md_render_mol_debug_kernel_info(md_render_mol_t r, uint32_t i, const char** name, md_gpu_kernel_info_t* info);
void md_render_debug_occlusion(float margin_radii, float min_px);

typedef struct pass_prof_t {
    md_gpu_stream_t s;
    md_tick_t       t;
    uint32_t        n;
    const char*     name[24];
    double          ms[24];
} pass_prof_t;

/* Each pass synced and timed on the host: the sync's own cost is in every number */
static void pass_hook(void* user, const char* pass) {
    pass_prof_t* p = user;
    md_gpu_stream_sync(p->s);
    md_tick_t now = md_tick_now();
    if (p->n < 24) { p->name[p->n] = pass; p->ms[p->n] = md_tick_to_milliseconds(now - p->t); p->n++; }
    p->t = now;
}

typedef struct timing_t { double best, median; } timing_t;

/* reps draws of the frame, each after `bump` (a content counter index, or -1 for none), each synced */
static timing_t time_draws(md_gpu_stream_t s, md_render_view_t v, md_render_frame_t* f, int bump, uint32_t reps, md_render_work_t* work) {
    double* t = malloc(reps * sizeof(double));
    md_render_status_t st = {0};
    for (uint32_t r = 0; r < reps; ++r) {
        if (bump >= 0) f->content[bump]++;
        md_tick_t t0 = md_tick_now();
        md_render_view_draw(s, v, f, &st);
        md_gpu_stream_sync(s);
        t[r] = md_tick_to_milliseconds(md_tick_now() - t0);
    }
    if (work) *work = st.work;
    qsort(t, reps, sizeof(double), cmp_d);
    timing_t out = { t[0], t[reps / 2] };
    free(t);
    return out;
}

static void print_times(const char* name, double* t, uint32_t reps) {
    qsort(t, reps, sizeof(double), cmp_d);
    printf("%s: %6.2f / %6.2f ms\n", name, t[0], t[reps / 2]);
}

/* The frame pipeline: a frame quantized on the host, uploaded into a slot of
   four (with and without the regroup, which reorders the three others), and
   the copy interpolated from 1, 2 or 4 resident frames. The box is periodic,
   so the minimum image is part of every interpolation. The frame is
   quantized straight into GPU-visible memory and read from there (on
   unified memory, nothing is copied); the upload from host memory, through
   the stream, is timed beside it. */
/* The best of `reps` frame updates from the given slots (velocities when dt) */
static double time_frames(md_gpu_stream_t s, md_render_mol_t mol, scene_t* sc, md_render_mol_atoms_t* atoms,
                          uint32_t count, const uint32_t slot[4], float dt, double* t, uint32_t reps) {
    atoms->frame_count = count;
    for (uint32_t k = 0; k < 4; ++k) atoms->frame_slot[k] = slot[k];
    for (uint32_t rr = 0; rr < reps + 2; ++rr) {
        md_render_mol_frame_weights(atoms->frame_w, atoms->frame_dw, count, (rr % 8) / 8.0f, dt, 0.5f);
        md_tick_t t0 = md_tick_now();
        md_render_mol_update(mol, s, &sc->topo, sc->dyn, atoms, 0, NULL);
        md_gpu_stream_sync(s);
        if (rr >= 2) t[rr - 2] = md_tick_to_milliseconds(md_tick_now() - t0);
    }
    qsort(t, reps, sizeof(double), cmp_d);
    return t[0];
}

static void bench_frames(md_gpu_stream_t s, md_render_mol_t mol, scene_t* sc, uint32_t reps) {
    const uint32_t n = sc->atoms;
    float cell[9] = { sc->box, 0, 0, 0, sc->box, 0, 0, 0, sc->box };
    md_gpu_mem_t up = md_gpu_malloc(s, MD_GPU_MEM_HOST_WRITE, md_render_mol_frame_size(&sc->topo));
    uint16_t* host_q = malloc((size_t)n * 6);
    const uint64_t size = md_render_mol_frames_size(&sc->topo, 4);
    md_gpu_addr_t base = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, size).gpu;
    md_render_mol_frames_t fr;
    double* t = malloc((reps + 2) * sizeof(double));
    if (!up.cpu || !base || !host_q || !md_render_mol_frames_init(&fr, &sc->topo, base, 4)) {
        printf("  frames: no memory for %.1f MB\n", size / 1048576.0);
        if (up.gpu) md_gpu_free(s, up.gpu);
        if (base) md_gpu_free(s, base);
        free(host_q); free(t);
        return;
    }
    md_render_mol_frame_t f, fh;
    md_bitset_t rset = sc->partial ? md_bitmask_view(&sc->resident) : md_bitset_range(0, n);
    const md_bitset_t* res = sc->partial ? &rset : NULL;
    for (uint32_t rr = 0; rr < reps; ++rr) {
        md_tick_t t0 = md_tick_now();
        md_render_mol_frame_quantize(&f, up.cpu, sc->host_pos, 0, n, cell, res);
        t[rr] = md_tick_to_milliseconds(md_tick_now() - t0);
    }
    f.q_gpu = up.gpu;
    print_times("frame quantize (host, one thread)   ", t, reps);
    if (sc->partial) {
        /* the read-ahead of a resident set's quantize, by its distance in
           windows of 64 atoms (0: none); the default is 16 */
        const uint32_t ahead[] = { 0, 2, 4, 8, 16, 32, 64 };
        printf("  quantize by read-ahead (windows of 64 atoms):");
        for (size_t ai = 0; ai < sizeof(ahead) / sizeof(ahead[0]); ++ai) {
            md_render_mol_debug_prefetch(ahead[ai]);
            for (uint32_t rr = 0; rr < reps; ++rr) {
                md_tick_t t0 = md_tick_now();
                md_render_mol_frame_quantize(&f, up.cpu, sc->host_pos, 0, n, cell, res);
                t[rr] = md_tick_to_milliseconds(md_tick_now() - t0);
            }
            qsort(t, reps, sizeof(double), cmp_d);
            printf(" %u: %.2f", ahead[ai], t[0]);
        }
        printf(" ms\n");
        md_render_mol_debug_prefetch(16);
        f.q_gpu = up.gpu;
    }
    md_render_mol_frame_quantize(&fh, host_q, sc->host_pos, 0, n, cell, res);

    for (uint32_t k = 0; k < 4; ++k) md_render_mol_frame_upload(mol, s, &sc->topo, 0, &fr, k, &f, 0);
    md_gpu_stream_sync(s);
    const struct { const char* name; const md_render_mol_frame_t* f; md_render_mol_update_flags_t flags; } ups[] = {
        { "frame upload (reorder)               ", &f,  0 },
        { "frame upload from host (copy, reorder)", &fh, 0 },
        { "frame upload (regroup, 3 permutes)   ", &f,  MD_RENDER_MOL_UPDATE_REGROUP },
    };
    for (int ui = 0; ui < 3; ++ui) {
        for (uint32_t rr = 0; rr < reps + 2; ++rr) {
            md_tick_t t0 = md_tick_now();
            md_render_mol_frame_upload(mol, s, &sc->topo, sc->dyn, &fr, rr % 4, ups[ui].f, ups[ui].flags);
            md_gpu_stream_sync(s);
            if (rr >= 2) t[rr - 2] = md_tick_to_milliseconds(md_tick_now() - t0);
        }
        print_times(ups[ui].name, t, reps);
    }
    pass_prof_t prof = { .s = s };
    md_gpu_stream_sync(s);
    prof.t = md_tick_now();
    md_render_mol_frame_upload_hooked(mol, s, &sc->topo, sc->dyn, &fr, 0, &f, MD_RENDER_MOL_UPDATE_REGROUP, pass_hook, &prof);
    printf("  passes, each synced (the first: the parameters and the clear):");
    for (uint32_t k = 0; k < prof.n; ++k) printf(" %s %.2f", prof.name[k], prof.ms[k]);
    printf("\n");

    const struct { const char* name; uint32_t count; float dt; } ints[] = {
        { "update from frames: nearest          ", 1, 0.0f },
        { "update from frames: linear           ", 2, 0.0f },
        { "update from frames: linear + velocity", 2, 1.0f },
        { "update from frames: cubic + velocity ", 4, 1.0f },
        { "update from frames: linear, again    ", 2, 0.0f },      /* the order is not the cause */
    };
    md_render_mol_atoms_t atoms = {0};
    atoms.frames = &fr;
    for (size_t ii = 0; ii < sizeof(ints) / sizeof(ints[0]); ++ii) {
        atoms.frame_count = ints[ii].count;
        for (uint32_t k = 0; k < 4; ++k) atoms.frame_slot[k] = k;
        for (uint32_t rr = 0; rr < reps + 2; ++rr) {
            md_render_mol_frame_weights(atoms.frame_w, atoms.frame_dw, ints[ii].count, (rr % 8) / 8.0f, ints[ii].dt, 0.5f);
            md_tick_t t0 = md_tick_now();
            md_render_mol_update(mol, s, &sc->topo, sc->dyn, &atoms, 0, NULL);
            md_gpu_stream_sync(s);
            if (rr >= 2) t[rr - 2] = md_tick_to_milliseconds(md_tick_now() - t0);
        }
        print_times(ints[ii].name, t, reps);
    }
    {
        /* the gather from fewer groups, more points a thread: whatever a
           thread does once (its arguments and the frames' parameters, read
           through pointers) spread over more points */
        const uint32_t s0123[4] = { 0, 1, 2, 3 };
        const uint32_t groups[] = { 65535, 8192, 2048, 640 };
        printf("  by groups (points a thread):    nearest  linear  cubic+vel\n");
        for (size_t gi = 0; gi < sizeof(groups) / sizeof(groups[0]); ++gi) {
            const uint32_t np = sc->topo.point_count;
            uint32_t g = groups[gi], full = (np + 63) / 64;
            if (g > full) g = full;
            md_render_mol_debug_frame_groups(g);
            const double a1 = time_frames(s, mol, sc, &atoms, 1, s0123, 0.0f, t, reps);
            const double a2 = time_frames(s, mol, sc, &atoms, 2, s0123, 0.0f, t, reps);
            const double a4 = time_frames(s, mol, sc, &atoms, 4, s0123, 1.0f, t, reps);
            printf("    %5u (%6.1f)               %6.2f  %6.2f  %6.2f\n", g, (double)np / (64.0 * g), a1, a2, a4);
        }
        md_render_mol_debug_frame_groups(0);
    }
    {
        /* the copy from positions again, beside the frames' (the same state of the GPU's clocks) */
        md_render_mol_atoms_t ap = {0};
        ap.position = sc->pos;
        for (uint32_t rr = 0; rr < reps + 2; ++rr) {
            md_tick_t t0 = md_tick_now();
            md_render_mol_update(mol, s, &sc->topo, sc->dyn, &ap, 0, NULL);
            md_gpu_stream_sync(s);
            if (rr >= 2) t[rr - 2] = md_tick_to_milliseconds(md_tick_now() - t0);
        }
        print_times("update from positions (the copy)     ", t, reps);
    }
    printf("  frames block: %.1f MB for 4 slots (%.1f bytes per resident atom), the upload buffer %.1f MB\n",
           size / 1048576.0, (double)size / (sc->topo.point_count ? sc->topo.point_count : 1), md_render_mol_frame_size(&sc->topo) / 1048576.0);
    md_gpu_free(s, base);
    md_gpu_free(s, up.gpu);
    free(host_q); free(t);
}

/* What reading the radii from the types' table costs, spacefill near and
   far, all options: the table, and a float radius per atom in the points
   instead (md_render_mol's override; the copy writes it). Three rounds of
   the two, their draws pooled: a GPU clock that drifts during the run moves
   both. */
static void bench_radii(md_gpu_stream_t s, md_render_t r, md_render_mol_t mol, scene_t* sc, const md_render_mol_atoms_t* atoms, uint32_t W, uint32_t H, uint32_t reps) {
    float* hr = malloc((size_t)sc->atoms * 4);
    md_gpu_addr_t arad = hr ? md_gpu_malloc(s, MD_GPU_MEM_DEVICE, (size_t)sc->atoms * 4).gpu : 0;
    if (!arad) { free(hr); return; }
    for (uint32_t i = 0; i < sc->atoms; ++i) hr[i] = (i % 3) ? 1.1f : 1.52f;   /* the scene's types' radii */
    md_gpu_upload(s, arad, hr, (size_t)sc->atoms * 4);
    md_gpu_stream_sync(s);
    free(hr);
    enum { ROUNDS = 3, VARS = 2, MEAS = 4 };
    const char* var[VARS] = { "table by type", "radius per atom" };
    const char* meas[MEAS] = { "raster (refit)", "raster (camera)", "state change", "refit pass" };
    const struct { const char* name; float dist; } looks[] = { { "near", 0.6f }, { "far", 6.0f } };
    double* t = malloc((size_t)VARS * MEAS * ROUNDS * reps * sizeof(double));
    printf("-- radii, spacefill (all options), %u rounds of %u draws: ms (best / median)\n", (unsigned)ROUNDS, reps);
    printf("  %-5s %-23s |", "", "");
    for (int m = 0; m < MEAS; ++m) printf(" %17s", meas[m]);
    printf("\n");
    for (size_t li = 0; li < sizeof(looks) / sizeof(looks[0]) && t; ++li) {
        md_render_view_t views[VARS];
        for (int vi = 0; vi < VARS; ++vi) views[vi] = md_render_view_create(r, "radii");
        for (uint32_t round = 0; round < ROUNDS; ++round) {
            for (int vi = 0; vi < VARS; ++vi) {
                md_render_mol_atoms_t at = *atoms;
                if (vi == 1) at.radius = arad;
                md_render_mol_update(mol, s, &sc->topo, sc->dyn, &at, 0, NULL);     /* the copy: w a type, or a radius */
                md_render_mol_rep_t rep = {0};
                rep.key = 1;
                rep.type = MD_RENDER_MOL_REP_SPACEFILL;
                rep.params.spacefill.radius_scale = 1.0f;
                rep.color.mode = MD_RENDER_COLOR_POINTS;
                rep.atom_color = sc->color;
                md_render_batch_t b[8];
                uint32_t nb = md_render_mol_emit(&sc->topo, sc->dyn, &at, &rep, b, 8);
                md_render_frame_t f;
                memset(&f, 0, sizeof(f));
                f.width = W; f.height = H;
                f.camera = camera_at(sc, looks[li].dist, (float)W / H);
                f.batches = b; f.batch_count = nb;
                f.jitter = MD_RENDER_JITTER_NONE; f.samples = 1;
                md_render_view_t v = views[vi];
                time_draws(s, v, &f, MD_RENDER_CONTENT_GEOMETRY, 2, NULL);          /* the new copy: a refit */
                double* o = t + ((size_t)vi * MEAS) * ROUNDS * reps + (size_t)round * reps;
                const size_t stride = (size_t)ROUNDS * reps;
                for (uint32_t k = 0; k < reps; ++k) {
                    f.content[MD_RENDER_CONTENT_GEOMETRY]++;
                    md_tick_t t0 = md_tick_now();
                    md_render_view_draw(s, v, &f, NULL);
                    md_gpu_stream_sync(s);
                    o[k] = md_tick_to_milliseconds(md_tick_now() - t0);
                }
                for (uint32_t k = 0; k < reps; ++k) {
                    f.camera.view[12] += 0.01f * sc->box / (float)reps * ((k & 1) ? -1.0f : 1.0f);
                    md_tick_t t0 = md_tick_now();
                    md_render_view_draw(s, v, &f, NULL);
                    md_gpu_stream_sync(s);
                    o[stride + k] = md_tick_to_milliseconds(md_tick_now() - t0);
                }
                for (uint32_t k = 0; k < reps; ++k) {
                    f.content[MD_RENDER_CONTENT_STATE]++;
                    md_tick_t t0 = md_tick_now();
                    md_render_view_draw(s, v, &f, NULL);
                    md_gpu_stream_sync(s);
                    o[2 * stride + k] = md_tick_to_milliseconds(md_tick_now() - t0);
                }
                for (uint32_t k = 0; k < reps; ++k) {
                    pass_prof_t prof = { .s = s };
                    f.content[MD_RENDER_CONTENT_GEOMETRY]++;
                    md_gpu_stream_sync(s);
                    prof.t = md_tick_now();
                    md_render_debug_pass_hook(pass_hook, &prof);
                    md_render_view_draw(s, v, &f, NULL);
                    md_render_debug_pass_hook(NULL, NULL);
                    md_gpu_stream_sync(s);
                    o[3 * stride + k] = 0.0;
                    for (uint32_t q = 0; q < prof.n; ++q) if (!strcmp(prof.name[q], "refit")) o[3 * stride + k] = prof.ms[q];
                }
            }
        }
        for (int vi = 0; vi < VARS; ++vi) {
            printf("  %-5s %-23s |", vi == 0 ? looks[li].name : "", var[vi]);
            for (int m = 0; m < MEAS; ++m) {
                double* o = t + ((size_t)vi * MEAS + m) * ROUNDS * reps;
                const size_t n = (size_t)ROUNDS * reps;
                qsort(o, n, sizeof(double), cmp_d);
                printf("   %6.2f / %6.2f", o[0], o[n / 2]);
            }
            printf("\n");
        }
        for (int vi = 0; vi < VARS; ++vi) md_render_view_destroy(s, views[vi]);
    }
    free(t);
    md_render_mol_update(mol, s, &sc->topo, sc->dyn, atoms, 0, NULL);
    md_gpu_stream_sync(s);
    md_gpu_free(s, arad);
}

static int run_gpu(int argc, char** argv) {
    uint32_t atoms_list[8] = { 100000, 1000000, 10000000 };
    uint32_t n_atoms = 3;
    uint32_t W = 1920, H = 1080, reps = 10;
    float smalls[8] = { 6 };
    uint32_t n_small = 1;
    cluster_src_t src = SRC_MOL;
    const char* only_look = NULL;     /* --look NAME: the looks whose name contains NAME */
    float resident = 1.0f;            /* --resident F: md_render_mol plans that fraction of the molecules */
    bool  resident_scattered = false; /* --resident-random F: ... spread at random */
    bool  cell_order = false;         /* --mol-order cells: the molecules in grid-cell order, not at random */
    float jitter = 0.0f;              /* --jitter A: with --mol-order cells, then moved by up to A along each axis */
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i], "--atoms") && i + 1 < argc) {
            n_atoms = 0;
            for (char* p = argv[++i]; *p && n_atoms < 8; ) { atoms_list[n_atoms++] = (uint32_t)strtoul(p, &p, 10); if (*p == ',') ++p; }
        } else if (!strcmp(argv[i], "--res") && i + 1 < argc) {
            sscanf(argv[++i], "%ux%u", &W, &H);
        } else if (!strcmp(argv[i], "--reps") && i + 1 < argc) {
            reps = (uint32_t)atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--implicit")) {
            src = SRC_IMPLICIT;
        } else if (!strcmp(argv[i], "--host-grid")) {
            src = SRC_HOST_GRID;
        } else if (!strcmp(argv[i], "--look") && i + 1 < argc) {
            only_look = argv[++i];
        } else if (!strcmp(argv[i], "--resident") && i + 1 < argc) {
            resident = strtof(argv[++i], NULL);
        } else if (!strcmp(argv[i], "--jitter") && i + 1 < argc) {
            jitter = (float)atof(argv[++i]);
        } else if (!strcmp(argv[i], "--mol-order") && i + 1 < argc) {
            cell_order = !strcmp(argv[++i], "cells");
        } else if (!strcmp(argv[i], "--resident-random") && i + 1 < argc) {
            resident = strtof(argv[++i], NULL);
            resident_scattered = true;
        } else if (!strcmp(argv[i], "--small") && i + 1 < argc) {
            n_small = 0;
            for (char* p = argv[++i]; *p && n_small < 8; ) { smalls[n_small++] = strtof(p, &p); if (*p == ',') ++p; }
        }
    }
    if (reps < 1) reps = 1;

    md_gpu_device_desc_t dd = {0};
    dd.label = "md_render_bench";
    md_gpu_device_t dev = md_gpu_device_create(&dd);
    if (!dev) { fprintf(stderr, "no GPU device: %s\n", md_gpu_last_error()); return 1; }
    md_gpu_device_info_t info;
    md_gpu_device_info(dev, &info);
    md_gpu_stream_t s = md_gpu_stream_default(dev, MD_GPU_STREAM_GRAPHICS);
    md_render_t r = md_render_create(dev, NULL);
    md_render_mol_t mol = md_render_mol_create(dev, NULL);
    if (!s || !r || !mol) { fprintf(stderr, "md_render could not start\n"); return 1; }

    printf("mdlib %s\n", MD_BENCH_SOURCE_REV);
    printf("device: %s (%s), 64-bit atomics: %s\n", info.name, info.driver, info.supports_atomic_int64 ? "yes" : "no (no compute raster)");
    printf("target %ux%u, %u draws per measurement, wall time per draw including the sync, ms (best / median)\n", W, H, reps);
    /* what the driver made of the kernels (registers, ...), where it says */
    {
        const char* kn; md_gpu_kernel_info_t ki;
        for (uint32_t i = 0; md_render_debug_kernel_info(r, i, &kn, &ki); ++i)
            if (ki.stats[0]) printf("  kernel md_render %-14s %s\n", kn, ki.stats);
        for (uint32_t i = 0; md_render_mol_debug_kernel_info(mol, i, &kn, &ki); ++i)
            if (ki.stats[0]) printf("  kernel md_render_mol %-10s %s\n", kn, ki.stats);
    }
    printf("\n");

    /* standin: licorice drawn where a stick is at least a pixel wide (2.5
       px/A for 0.2 A), spheres at 0.4 of the atoms' radii beyond (md_render_mol
       only: two representations in adjacent size windows) */
    const struct { const char* name; bool licorice; float dist; float standin; } looks[] = {
        { "spacefill, near", false, 0.6f, 0.0f },
        { "spacefill, far",  false, 6.0f, 0.0f },
        { "licorice, near",  true,  0.6f, 0.0f },
        { "licorice from 2.5 px/A, spheres beyond, near", true, 0.6f, 2.5f },
    };
    /* md_render's occlusion defaults: margin 2 radii, elements up to 2 px
       untested. Rows marked "x:" would be experiments with other settings
       (md_render_debug_occlusion); none at the moment. */
    const struct { const char* name; uint32_t options; float margin, min_px; } opts[] = {
        { "all",                     0, 2.0f, 2.0f },
        { "no occlusion",            MD_RENDER_NO_OCCLUSION_CULLING, 2.0f, 2.0f },
        { "no compute raster",       MD_RENDER_NO_COMPUTE_RASTER, 2.0f, 2.0f },
        { "neither (plain raster)",  MD_RENDER_NO_OCCLUSION_CULLING | MD_RENDER_NO_COMPUTE_RASTER, 2.0f, 2.0f },
    };

    for (uint32_t ai = 0; ai < n_atoms; ++ai) {
        scene_t sc;
        if (!scene_make(s, &sc, atoms_list[ai] / 3, src, resident, resident_scattered, cell_order, jitter)) { fprintf(stderr, "scene of %u atoms failed (memory?)\n", atoms_list[ai]); continue; }
        printf("== %u atoms, %u bonds, box %.0f A, %s\n", sc.atoms, sc.bonds, sc.box,
               src == SRC_IMPLICIT ? "implicit clusters (runs of 64 in grid-cell order)" :
               src == SRC_HOST_GRID ? "clusters by grid cell, built on the host" : (cell_order ? "md_render_mol (molecules in grid-cell order, GPU grid)" : "md_render_mol (molecules in random order, GPU grid)"));
        if (cell_order && jitter > 0.0f) printf("   then each molecule moved by up to %.1f A along each axis\n", jitter);
        if (src == SRC_MOL) printf("   %u atoms resident (%.1f %%%s), %u bonds; md_render_mol_plan %.1f ms (host, one thread), write %.1f ms\n",
                                   sc.topo.point_count, 100.0 * sc.topo.point_count / sc.atoms,
                                   !sc.partial ? "" : resident_scattered ? ", molecules at random" : ", one block of molecules",
                                   sc.topo.bond_count, sc.plan_ms, sc.write_ms);
        md_render_mol_atoms_t atoms = {0};
        atoms.position = sc.pos; atoms.selection = sc.sel;
        if (src == SRC_MOL) {
            /* what a position change costs on top of the raster: the copy, and the copy after a regroup */
            const struct { const char* name; md_render_mol_update_flags_t flags; } ups[] = {
                { "md_render_mol_update (copy)          ", 0 },
                { "md_render_mol_update (regroup + copy)", MD_RENDER_MOL_UPDATE_REGROUP },
            };
            double* tu = malloc(reps * sizeof(double));
            for (int ui = 1; ui >= 0; --ui) {          /* the regroup first: the copy then times the grouping it leaves */
                for (uint32_t rr = 0; rr < reps + 2; ++rr) {
                    md_tick_t t0 = md_tick_now();
                    md_render_mol_update(mol, s, &sc.topo, sc.dyn, &atoms, ups[ui].flags, NULL);
                    md_gpu_stream_sync(s);
                    if (rr >= 2) tu[rr - 2] = md_tick_to_milliseconds(md_tick_now() - t0);
                }
                qsort(tu, reps, sizeof(double), cmp_d);
                printf("%s: %6.2f / %6.2f ms\n", ups[ui].name, tu[0], tu[reps / 2]);
            }
            if (sc.host_pos && sc.plan.data) {
                /* the regroup's other cases: after the molecules moved (by up to
                   12 A along each axis since the last grouping: positions that
                   alternate between the scene's and a copy so moved), and the
                   first after the plan (from molecule order: scattered once;
                   the topology written again before each) */
                const float J = 12.0f;
                float* hj = malloc((size_t)sc.atoms * 12);
                uint32_t js = 4242;
                #define JRND() (((js = js * 1664525u + 1013904223u) >> 8) / 16777216.0f)
                for (uint32_t m = 0; m < sc.atoms / 3; ++m) {
                    const float d[3] = { (2.0f * JRND() - 1.0f) * J, (2.0f * JRND() - 1.0f) * J, (2.0f * JRND() - 1.0f) * J };
                    for (int k = 0; k < 9; ++k) hj[9 * (size_t)m + k] = sc.host_pos[9 * (size_t)m + k] + d[k % 3];
                }
                #undef JRND
                md_gpu_addr_t pj = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, (size_t)sc.atoms * 12).gpu;
                md_gpu_upload(s, pj, hj, (size_t)sc.atoms * 12);
                free(hj);
                md_render_mol_atoms_t aj = atoms;
                aj.position = pj;
                for (uint32_t rr = 0; rr < reps + 2; ++rr) {
                    md_tick_t t0 = md_tick_now();
                    md_render_mol_update(mol, s, &sc.topo, sc.dyn, (rr & 1) ? &aj : &atoms, MD_RENDER_MOL_UPDATE_REGROUP, NULL);
                    md_gpu_stream_sync(s);
                    if (rr >= 2) tu[rr - 2] = md_tick_to_milliseconds(md_tick_now() - t0);
                }
                qsort(tu, reps, sizeof(double), cmp_d);
                printf("md_render_mol_update (regroup + copy, molecules moved up to %.0f A): %6.2f / %6.2f ms\n", J, tu[0], tu[reps / 2]);
                const uint32_t nf = reps < 3 ? reps : 3;
                for (uint32_t rr = 0; rr < nf; ++rr) {
                    md_render_mol_write(s, &sc.topo, &sc.plan, sc.topo.base, NULL);
                    md_gpu_stream_sync(s);
                    md_tick_t t0 = md_tick_now();
                    md_render_mol_update(mol, s, &sc.topo, sc.dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, NULL);
                    md_gpu_stream_sync(s);
                    tu[rr] = md_tick_to_milliseconds(md_tick_now() - t0);
                }
                qsort(tu, nf, sizeof(double), cmp_d);
                printf("md_render_mol_update (the first regroup after the plan + copy): %6.2f / %6.2f ms\n", tu[0], tu[nf / 2]);
                md_render_mol_update(mol, s, &sc.topo, sc.dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, NULL);
                md_gpu_stream_sync(s);
                md_gpu_free(s, pj);
            }
            free(tu);
            /* where the regroup's time goes: each pass synced (so each carries a sync) */
            pass_prof_t prof = { .s = s };
            md_gpu_stream_sync(s);
            prof.t = md_tick_now();
            md_render_mol_update_hooked(mol, s, &sc.topo, sc.dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, pass_hook, &prof);
            printf("  passes, each synced:");
            for (uint32_t k = 0; k < prof.n; ++k) printf(" %s %.2f", prof.name[k], prof.ms[k]);
            {
                /* the write pass leaving out the atoms' order, the bonds, then
                   both (the records it always writes: the next regroup reads
                   them); best of 3 each. The topology is wrong meanwhile: a full
                   regroup restores it. */
                const struct { uint32_t mask; const char* name; } wp[] = {
                    { 0, "all" }, { 1, "no atom_order" }, { 4, "no bonds" }, { 5, "records only" },
                };
                printf("\n  the write pass by what it writes (ms):");
                for (size_t wi = 0; wi < sizeof(wp) / sizeof(wp[0]); ++wi) {
                    md_render_mol_debug_write(wp[wi].mask);
                    double best = 1e30;
                    for (int rep = 0; rep < 3; ++rep) {
                        pass_prof_t pw = { .s = s };
                        md_gpu_stream_sync(s);
                        pw.t = md_tick_now();
                        md_render_mol_update_hooked(mol, s, &sc.topo, sc.dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, pass_hook, &pw);
                        for (uint32_t k = 0; k < pw.n; ++k) if (!strcmp(pw.name[k], "write") && pw.ms[k] < best) best = pw.ms[k];
                    }
                    printf("%s %s %.2f", wi ? "," : "", wp[wi].name, best);
                }
                md_render_mol_debug_write(0);
                md_render_mol_update(mol, s, &sc.topo, sc.dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, NULL);
                md_gpu_stream_sync(s);
            }
            printf("\n  %u grid cells (capacity), %.1f MB topology, %.1f MB per-frame block (%.1f MB with velocities)\n",
                   sc.topo.grid_cells, sc.topo.size / 1048576.0, sc.topo.dynamic_size / 1048576.0, sc.topo.dynamic_size_velocity / 1048576.0);
            bench_frames(s, mol, &sc, reps);
            /* the copy drawn below: from the positions, as before */
            md_render_mol_update(mol, s, &sc.topo, sc.dyn, &atoms, 0, NULL);
            md_gpu_stream_sync(s);
        }
        for (size_t li = 0; li < sizeof(looks) / sizeof(looks[0]); ++li) {
            if (only_look && !strstr(looks[li].name, only_look)) continue;
            if (looks[li].standin > 0.0f && src != SRC_MOL) continue;
            md_render_batch_t b[8];
            memset(b, 0, sizeof(b));
            uint32_t nb = 0;
            if (src == SRC_MOL) {
                md_render_mol_rep_t rep = {0};
                rep.key = 1;
                rep.color.mode = MD_RENDER_COLOR_POINTS;
                rep.atom_color = sc.color;
                if (looks[li].licorice) { rep.type = MD_RENDER_MOL_REP_LICORICE; rep.params.licorice.radius = 0.2f; }
                else                    { rep.type = MD_RENDER_MOL_REP_SPACEFILL; rep.params.spacefill.radius_scale = 1.0f; }
                rep.lod_min = looks[li].standin;
                nb = md_render_mol_emit(&sc.topo, sc.dyn, &atoms, &rep, b, 8);
                if (looks[li].standin > 0.0f) {
                    md_render_mol_rep_t far = rep;
                    far.key = 2;
                    far.type = MD_RENDER_MOL_REP_SPACEFILL;
                    far.params.spacefill.radius_scale = 0.4f;
                    far.lod_min = 0.0f;
                    far.lod_max = looks[li].standin;
                    nb += md_render_mol_emit(&sc.topo, sc.dyn, &atoms, &far, b + nb, 8 - nb);
                }
            }
            md_render_points_t pts = {0};
            pts.count = sc.atoms; pts.position = sc.pos; pts.radius = sc.rad; pts.color = sc.color; pts.state[0] = sc.sel;
            if (src != SRC_MOL) {
                b[nb].key = 1; b[nb].kind = MD_RENDER_SPHERES; b[nb].count = sc.atoms; b[nb].points = pts;
                if (src == SRC_HOST_GRID) { b[nb].cluster_offset = sc.atom_clusters; b[nb].cluster_count = sc.clusters; }
                if (looks[li].licorice) b[nb].style.spheres.radius_add = 0.2f; else b[nb].style.spheres.radius_scale = 1.0f;
                nb++;
                if (looks[li].licorice) {
                    b[nb].key = 2; b[nb].kind = MD_RENDER_CAPSULES; b[nb].count = sc.bonds; b[nb].points = pts;
                    b[nb].index = sc.bond_idx; b[nb].style.capsules.radius = 0.2f;
                    if (src == SRC_HOST_GRID) { b[nb].cluster_offset = sc.bond_clusters; b[nb].cluster_count = sc.clusters; }
                    nb++;
                }
            }
            printf("-- %s\n", looks[li].name);
            printf("  %-24s %6s | %17s %17s %17s %17s | %9s %9s %9s %9s %9s %9s | %8s %8s %8s\n", "options", "small", "raster (refit)", "raster (camera)", "resolve", "mask",
                   "outside", "lod", "occluded", "hw", "compute", "retested", "clusters", "expanded", "c.occl");
            for (size_t oi = 0; oi < sizeof(opts) / sizeof(opts[0]); ++oi) {
                for (uint32_t pi = 0; pi < (opts[oi].options & MD_RENDER_NO_COMPUTE_RASTER ? 1u : n_small); ++pi) {
                    md_render_debug_occlusion(opts[oi].margin, opts[oi].min_px);
                    md_render_view_t v = md_render_view_create(r, "bench");
                    md_render_frame_t f;
                    memset(&f, 0, sizeof(f));
                    f.width = W; f.height = H;
                    f.camera = camera_at(&sc, looks[li].dist, (float)W / H);
                    f.batches = b; f.batch_count = nb;
                    f.jitter = MD_RENDER_JITTER_NONE; f.samples = 1;
                    f.options = opts[oi].options;
                    f.small_pixels = smalls[pi];
                    /* warm up: the first raster has no pyramid */
                    time_draws(s, v, &f, MD_RENDER_CONTENT_GEOMETRY, 2, NULL);
                    timing_t tr = time_draws(s, v, &f, MD_RENDER_CONTENT_GEOMETRY, reps, NULL);
                    md_render_stats_t st = {0};
                    md_render_view_stats(v, &st);
                    /* the camera moving: a raster without a refit (one camera step per draw) */
                    double* tc = malloc(reps * sizeof(double));
                    for (uint32_t rr = 0; rr < reps; ++rr) {
                        f.camera.view[12] += 0.01f * sc.box / (float)reps * ((rr & 1) ? -1.0f : 1.0f);
                        md_tick_t t0 = md_tick_now();
                        md_render_view_draw(s, v, &f, NULL);
                        md_gpu_stream_sync(s);
                        tc[rr] = md_tick_to_milliseconds(md_tick_now() - t0);
                    }
                    qsort(tc, reps, sizeof(double), cmp_d);
                    timing_t tcam = { tc[0], tc[reps / 2] };
                    free(tc);
                    timing_t tv = time_draws(s, v, &f, MD_RENDER_CONTENT_APPEARANCE, reps, NULL);
                    timing_t tm = time_draws(s, v, &f, MD_RENDER_CONTENT_STATE, reps, NULL);
                    char small_s[16];
                    if (opts[oi].options & MD_RENDER_NO_COMPUTE_RASTER) snprintf(small_s, sizeof(small_s), "-");
                    else snprintf(small_s, sizeof(small_s), "%.0f", smalls[pi]);
                    printf("  %-24s %6s | %7.2f / %7.2f %7.2f / %7.2f %7.2f / %7.2f %7.2f / %7.2f | %9u %9u %9u %9u %9u %9u | %8u %8u %8u\n", opts[oi].name, small_s,
                           tr.best, tr.median, tcam.best, tcam.median, tv.best, tv.median, tm.best, tm.median,
                           st.outside, st.lod, st.occluded, st.hardware[0] + st.hardware[1], st.compute[0] + st.compute[1], st.occlusion_retested,
                           st.clusters, st.clusters_expanded[0] + st.clusters_expanded[1], st.clusters_occluded);
                    if (oi == 0 && pi == 0) {
                        /* where a camera raster's time goes, pass by pass (each synced) */
                        pass_prof_t prof = { .s = s };
                        f.camera.view[12] += 0.001f * sc.box;
                        md_gpu_stream_sync(s);
                        prof.t = md_tick_now();
                        md_render_debug_pass_hook(pass_hook, &prof);
                        md_render_view_draw(s, v, &f, NULL);
                        md_render_debug_pass_hook(NULL, NULL);
                        md_gpu_stream_sync(s);
                        printf("    camera raster, passes each synced:");
                        for (uint32_t k = 0; k < prof.n; ++k) printf(" %s %.2f", prof.name[k], prof.ms[k]);
                        md_render_stats_t ps = {0};
                        md_render_view_stats(v, &ps);
                        printf("\n    by phase: clusters expanded %u + %u, hardware %u + %u, compute %u + %u, occluded %u (retested %u)\n",
                               ps.clusters_expanded[0], ps.clusters_expanded[1], ps.hardware[0], ps.hardware[1], ps.compute[0], ps.compute[1],
                               ps.occluded, ps.occlusion_retested);
                        /* the same for a raster after a change of positions (the refit
                           and the state lists on top), and for a change of the state
                           bits alone (the selection: their refit, lists and the mask) */
                        const struct { const char* name; int bump; } more[] = {
                            { "refit raster", MD_RENDER_CONTENT_GEOMETRY },
                            { "state change", MD_RENDER_CONTENT_STATE },
                        };
                        for (int mi = 0; mi < 2; ++mi) {
                            prof.n = 0;
                            f.content[more[mi].bump]++;
                            md_gpu_stream_sync(s);
                            prof.t = md_tick_now();
                            md_render_debug_pass_hook(pass_hook, &prof);
                            md_render_view_draw(s, v, &f, NULL);
                            md_render_debug_pass_hook(NULL, NULL);
                            md_gpu_stream_sync(s);
                            printf("    %s, passes each synced:", more[mi].name);
                            for (uint32_t k = 0; k < prof.n; ++k) printf(" %s %.2f", prof.name[k], prof.ms[k]);
                            printf("\n");
                        }
                        if (getenv("MD_BENCH_SEQ")) {
                            float base12 = f.camera.view[12];
                            const char* step[] = { "c0 (geometry)", "c0 again (geometry)", "c1 = c0 + shift", "c1 again (geometry)", "c0 (back)", "c0 again (geometry)" };
                            const float sh[] = { 0, 0, 1, 1, 0, 0 };
                            const int  geo[] = { 1, 1, 0, 1, 0, 1 };
                            for (int q = 0; q < 6; ++q) {
                                f.camera.view[12] = base12 + sh[q] * (getenv("MD_BENCH_SHIFT") ? (float)atof(getenv("MD_BENCH_SHIFT")) : 0.001f) * sc.box;
                                if (geo[q]) f.content[MD_RENDER_CONTENT_GEOMETRY]++;
                                md_render_view_draw(s, v, &f, NULL);
                                md_gpu_stream_sync(s);
                                md_render_view_stats(v, &ps);
                                printf("    %-22s: clusters %u + %u, hw %u + %u, compute %u + %u, occluded %u (retested %u)\n", step[q],
                                       ps.clusters_expanded[0], ps.clusters_expanded[1], ps.hardware[0], ps.hardware[1], ps.compute[0], ps.compute[1],
                                       ps.occluded, ps.occlusion_retested);
                            }
                        }
                    }
                    md_render_view_destroy(s, v);
                    md_render_debug_occlusion(2.0f, 2.0f);
                }
            }
        }
        if (src == SRC_MOL && (!only_look || strstr("radii", only_look))) bench_radii(s, r, mol, &sc, &atoms, W, H, reps);
        md_gpu_memory_stats_t ms;
        if (md_gpu_memory_stats(dev, MD_GPU_MEM_DEVICE, &ms)) printf("device memory peak: %.0f MB\n", ms.bytes_peak_in_use / 1048576.0);
        printf("\n");
        scene_free(s, &sc);
        md_gpu_stream_sync(s);
    }
    md_render_mol_destroy(mol);
    md_render_destroy(r);
    md_gpu_device_destroy(dev);
    return 0;
}
#endif

int main(int argc, char** argv) {
    if (argc >= 2 && !strcmp(argv[1], "clusters")) return run_clusters(argc - 2, argv + 2);
#if MD_ENABLE_GPU
    if (argc >= 2 && !strcmp(argv[1], "gpu")) return run_gpu(argc - 2, argv + 2);
#endif
    fprintf(stderr, "usage:\n  md_render_bench clusters <structure> [<trajectory>] [--frames N] [--size N]\n"
#if MD_ENABLE_GPU
                    "  md_render_bench gpu [--atoms N[,N...]] [--res WxH] [--reps N] [--small P[,P...]] [--host-grid | --implicit] [--look NAME] [--resident F | --resident-random F] [--mol-order cells [--jitter A]]\n"
#endif
           );
    return 1;
}

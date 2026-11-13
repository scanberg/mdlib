// md_render_snap: a structure through md_render_mol and md_render, written as
// PNG images. Headless: no window, so it runs where md_render_bench does
// (lavapipe included), for looking at what the renderer draws on real data:
// representations, fades, masks.
//
//   md_render_snap <structure> [options]
//     --coords FILE        coordinates for it (a coordinate file, or a trajectory's first frame)
//     --rep TYPE[:P]       spacefill[:radius scale], licorice[:radius], ball[:ball scale],
//                          bonds[:radius], cartoon[:scale]; repeatable, in draw order
//                          (default spacefill). The cartoon takes the state's secondary
//                          structure (computed when it has none) and its own colours:
//                          helix red, sheet yellow, coil grey
//     --fade MODE S0 S1    the last --rep fades out across scales S0 > S1 (pixels per A),
//                          MODE thin, shrink or both
//     --fade-in MODE S0 S1 the last --rep fades in across S0 > S1
//     --view NAME          front, back, side, top or iso; repeatable (default front)
//     --zoom Z             camera distance, 1 = the scene just fits (default 1)
//     --replicate N        N x N x N copies along the cell (its box, or the bounds), as instances
//     --select A:B         atoms A .. B-1 selected (the mask)
//     --time N             N draws after a geometry change per view, best and median ms,
//                          then one with every pass synced and timed
//     --res WxH            (default 1280x720)
//     --out PREFIX         files PREFIX_<view>.png (default snap)
//
// The image: the resolved base colour shaded by its normal (a light from over
// the camera's shoulder), the selection's visible part tinted and its hidden
// part outlined. It is a look at the targets, not viamd's shading.

#include <md_system.h>
#include <md_system_state.h>
#include <md_infer.h>
#include <md_io.h>
#include <md_unitcell.h>
#include <md_element.h>
#include <md_backbone.h>
#include <md_render.h>
#include <md_render_mol.h>
#include <core/md_gpu.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>
#include <core/md_os.h>

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bench_rev.h"                    // MD_BENCH_SOURCE_REV

#if MD_ENABLE_GPU

/* ---- PNG: stored (uncompressed) deflate, no dependencies ---------------------- */

static uint32_t crc_table[256];
static void crc_init(void) {
    for (uint32_t n = 0; n < 256; ++n) {
        uint32_t c = n;
        for (int k = 0; k < 8; ++k) c = (c & 1) ? 0xEDB88320u ^ (c >> 1) : c >> 1;
        crc_table[n] = c;
    }
}
static uint32_t crc_update(uint32_t c, const uint8_t* p, size_t n) {
    for (size_t i = 0; i < n; ++i) c = crc_table[(c ^ p[i]) & 0xFF] ^ (c >> 8);
    return c;
}
static void put_be32(uint8_t* p, uint32_t v) { p[0] = (uint8_t)(v >> 24); p[1] = (uint8_t)(v >> 16); p[2] = (uint8_t)(v >> 8); p[3] = (uint8_t)v; }
static void png_chunk(FILE* f, const char* type, const uint8_t* data, uint32_t len) {
    uint8_t h[8];
    put_be32(h, len);
    memcpy(h + 4, type, 4);
    fwrite(h, 1, 8, f);
    if (len) fwrite(data, 1, len, f);
    uint32_t c = crc_update(0xFFFFFFFFu, (const uint8_t*)type, 4);
    c = crc_update(c, data, len) ^ 0xFFFFFFFFu;
    uint8_t t[4];
    put_be32(t, c);
    fwrite(t, 1, 4, f);
}
/* rgb: w * h * 3 bytes, rows top to bottom */
static bool write_png(const char* path, const uint8_t* rgb, uint32_t w, uint32_t h) {
    FILE* f = fopen(path, "wb");
    if (!f) return false;
    static const uint8_t sig[8] = { 137, 80, 78, 71, 13, 10, 26, 10 };
    fwrite(sig, 1, 8, f);
    uint8_t ihdr[13];
    put_be32(ihdr, w); put_be32(ihdr + 4, h);
    ihdr[8] = 8; ihdr[9] = 2; ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;      /* 8 bit RGB */
    png_chunk(f, "IHDR", ihdr, 13);
    const size_t row = (size_t)w * 3 + 1, raw = row * h;
    const size_t blocks = (raw + 65534) / 65535;
    const size_t zlen = 2 + raw + blocks * 5 + 4;
    uint8_t* z = malloc(zlen);
    uint8_t* src = malloc(raw);
    if (!z || !src) { free(z); free(src); fclose(f); return false; }
    for (uint32_t y = 0; y < h; ++y) { src[y * row] = 0; memcpy(src + y * row + 1, rgb + (size_t)y * w * 3, (size_t)w * 3); }
    size_t o = 0;
    z[o++] = 0x78; z[o++] = 0x01;
    uint32_t a = 1, b = 0;
    for (size_t i = 0; i < raw; ++i) { a = (a + src[i]) % 65521u; b = (b + a) % 65521u; }
    for (size_t pos = 0; pos < raw; pos += 65535) {
        uint32_t n = (uint32_t)(raw - pos < 65535 ? raw - pos : 65535);
        z[o++] = pos + n >= raw ? 1 : 0;
        z[o++] = (uint8_t)n; z[o++] = (uint8_t)(n >> 8);
        z[o++] = (uint8_t)~n; z[o++] = (uint8_t)(~n >> 8);
        memcpy(z + o, src + pos, n);
        o += n;
    }
    put_be32(z + o, (b << 16) | a);
    o += 4;
    png_chunk(f, "IDAT", z, (uint32_t)o);
    png_chunk(f, "IEND", NULL, 0);
    free(z); free(src);
    return fclose(f) == 0;
}

/* ---- Image ---------------------------------------------------------------------- */

static float half_to_float(uint16_t h) {
    uint32_t s = (uint32_t)(h & 0x8000u) << 16, e = (h >> 10) & 0x1Fu, m = h & 0x3FFu, x;
    if (e == 0) {
        if (m == 0) x = s;
        else { e = 127 - 15 + 1; while (!(m & 0x400u)) { m <<= 1; e--; } m &= 0x3FFu; x = s | (e << 23) | (m << 13); }
    } else if (e == 31) x = s | 0x7F800000u | (m << 13);
    else x = s | ((e + 127 - 15) << 23) | (m << 13);
    float f;
    memcpy(&f, &x, 4);
    return f;
}
static float srgb_to_lin(float c) { return c <= 0.04045f ? c / 12.92f : powf((c + 0.055f) / 1.055f, 2.4f); }
static float lin_to_srgb(float c) { c = c < 0 ? 0 : c > 1 ? 1 : c; return c <= 0.0031308f ? c * 12.92f : 1.055f * powf(c, 1.0f / 2.4f) - 0.055f; }

static bool read_tex(md_gpu_stream_t s, md_gpu_texture_t t, void* dst, size_t size) {
    md_gpu_tex_region_t r = {0};
    r.extent[2] = 1;
    if (md_gpu_texture_region_size(t, &r) != size) return false;
    md_gpu_mem_t rb = md_gpu_malloc(s, MD_GPU_MEM_HOST_READ, size);
    if (!rb.cpu) return false;
    bool ok = md_gpu_copy_from_texture(s, rb.gpu, t, &r);
    md_gpu_stream_sync(s);
    if (ok) memcpy(dst, rb.cpu, size);
    md_gpu_free(s, rb.gpu);
    return ok;
}

/* The shaded image of a view's targets */
static bool compose(md_gpu_stream_t s, md_render_view_t v, uint32_t W, uint32_t H, uint8_t* rgb) {
    md_render_targets_t t = md_render_view_targets(v);
    uint32_t* col = malloc((size_t)W * H * 4);
    uint32_t* nrm = malloc((size_t)W * H * 4);
    uint32_t* msk = malloc((size_t)W * H * 4);
    bool ok = col && nrm && msk && read_tex(s, t.color, col, (size_t)W * H * 4) && read_tex(s, t.normal, nrm, (size_t)W * H * 4) &&
              read_tex(s, t.mask, msk, (size_t)W * H * 4);
    if (ok) {
        float L[3] = { -0.35f, 0.45f, 0.82f };              /* view space, towards the viewer and up-left */
        float ll = sqrtf(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
        for (int k = 0; k < 3; ++k) L[k] /= ll;
        for (uint32_t y = 0; y < H; ++y) {
            for (uint32_t x = 0; x < W; ++x) {
                const size_t i = (size_t)y * W + x;
                uint8_t* o = rgb + i * 3;
                float out[3] = { 0.13f, 0.13f, 0.15f };     /* background, sRGB */
                if (col[i] >> 24) {
                    float ox = half_to_float((uint16_t)(nrm[i] & 0xFFFFu)), oy = half_to_float((uint16_t)(nrm[i] >> 16));
                    float n[3] = { ox, oy, 1.0f - fabsf(ox) - fabsf(oy) };
                    if (n[2] < 0.0f) {
                        float px = (1.0f - fabsf(oy)) * (ox >= 0.0f ? 1.0f : -1.0f);
                        float py = (1.0f - fabsf(ox)) * (oy >= 0.0f ? 1.0f : -1.0f);
                        n[0] = px; n[1] = py;
                    }
                    float nl = sqrtf(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
                    float d = nl > 0 ? (n[0] * L[0] + n[1] * L[1] + n[2] * L[2]) / nl : 0.0f;
                    float shade = 0.22f + 0.78f * (d > 0 ? d : 0.0f);
                    for (int c = 0; c < 3; ++c) out[c] = lin_to_srgb(srgb_to_lin(((col[i] >> (8 * c)) & 0xFF) / 255.0f) * shade);
                }
                /* the selection: its visible part tinted, its outline (hidden too) drawn */
                const uint32_t m = msk[i];
                const bool vis = (m & 0xFFu) != 0, any = ((m >> 8) & 0xFFu) != 0;
                if (vis) for (int c = 0; c < 3; ++c) out[c] = out[c] * 0.6f + (c == 2 ? 0.1f : 0.95f) * 0.4f;
                if (any) {
                    bool edge = false;
                    const int dx[4] = { 1, -1, 0, 0 }, dy[4] = { 0, 0, 1, -1 };
                    for (int k = 0; k < 4 && !edge; ++k) {
                        int xx = (int)x + dx[k], yy = (int)y + dy[k];
                        if (xx < 0 || yy < 0 || xx >= (int)W || yy >= (int)H) continue;
                        edge = ((msk[(size_t)yy * W + xx] >> 8) & 0xFFu) == 0;
                    }
                    if (edge) { out[0] = 1.0f; out[1] = 0.85f; out[2] = 0.1f; }
                }
                for (int c = 0; c < 3; ++c) { float q = out[c] * 255.0f + 0.5f; o[c] = (uint8_t)(q < 0 ? 0 : q > 255 ? 255 : q); }
            }
        }
    }
    free(col); free(nrm); free(msk);
    return ok;
}

/* md_render's per-pass hook: in md_render.c, not in its header (the bench's) */
void md_render_debug_pass_hook(void (*hook)(void* user, const char* pass), void* user);

typedef struct pass_prof_t {
    md_gpu_stream_t s;
    md_tick_t       t;
    uint32_t        n;
    const char*     name[32];
    double          ms[32];
} pass_prof_t;

/* Each pass synced and timed on the host: the sync's own cost is in every number */
static void pass_hook(void* user, const char* pass) {
    pass_prof_t* p = user;
    md_gpu_stream_sync(p->s);
    md_tick_t now = md_tick_now();
    if (p->n < 32) { p->name[p->n] = pass; p->ms[p->n] = md_tick_to_milliseconds(now - p->t); p->n++; }
    p->t = now;
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

/* ---- Scene ----------------------------------------------------------------------- */

typedef struct rep_arg_t {
    md_render_mol_rep_type_t type;
    float                    param;      /* < 0: the type's default */
    md_render_fade_t         fade;
} rep_arg_t;

static bool parse_fade_mode(const char* s, uint32_t* mode) {
    if (!strcmp(s, "thin"))   { *mode = MD_RENDER_FADE_THIN; return true; }
    if (!strcmp(s, "shrink")) { *mode = MD_RENDER_FADE_SHRINK; return true; }
    if (!strcmp(s, "both"))   { *mode = MD_RENDER_FADE_THIN | MD_RENDER_FADE_SHRINK; return true; }
    return false;
}

static void mat_ident(float m[16]) { memset(m, 0, 64); m[0] = m[5] = m[10] = m[15] = 1.0f; }

/* A camera looking at `ctr` from direction `dir` (unit, from the target to the
   eye) at `dist`, reverse Z, infinite far plane */
static md_render_camera_t camera_look(const double ctr[3], const double dir[3], double dist, float aspect, double near_z) {
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

static int run(int argc, char** argv) {
    if (argc < 1) return 2;
    const char* structure = argv[0];
    const char* coords = NULL;
    const char* out_prefix = "snap";
    uint32_t W = 1280, H = 720, replicate = 1, timed = 0;
    double zoom = 1.0;
    int64_t sel_beg = -1, sel_end = -1;
    rep_arg_t reps[8];
    uint32_t n_reps = 0;
    const char* views[8];
    uint32_t n_views = 0;
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (!strcmp(a, "--coords") && i + 1 < argc) coords = argv[++i];
        else if (!strcmp(a, "--out") && i + 1 < argc) out_prefix = argv[++i];
        else if (!strcmp(a, "--res") && i + 1 < argc) { if (sscanf(argv[++i], "%ux%u", &W, &H) != 2) return 2; }
        else if (!strcmp(a, "--zoom") && i + 1 < argc) zoom = atof(argv[++i]);
        else if (!strcmp(a, "--replicate") && i + 1 < argc) replicate = (uint32_t)atoi(argv[++i]);
        else if (!strcmp(a, "--time") && i + 1 < argc) timed = (uint32_t)atoi(argv[++i]);
        else if (!strcmp(a, "--select") && i + 1 < argc) { long long x, y; if (sscanf(argv[++i], "%lld:%lld", &x, &y) != 2) return 2; sel_beg = x; sel_end = y; }
        else if (!strcmp(a, "--view") && i + 1 < argc) { if (n_views < 8) views[n_views++] = argv[++i]; else ++i; }
        else if (!strcmp(a, "--rep") && i + 1 < argc) {
            if (n_reps == 8) return 2;
            char name[32] = {0};
            float p = -1.0f;
            const char* s = argv[++i];
            const char* colon = strchr(s, ':');
            size_t len = colon ? (size_t)(colon - s) : strlen(s);
            if (len >= sizeof(name)) return 2;
            memcpy(name, s, len);
            if (colon) p = (float)atof(colon + 1);
            rep_arg_t r;
            memset(&r, 0, sizeof(r));
            r.param = p;
            if      (!strcmp(name, "spacefill")) r.type = MD_RENDER_MOL_REP_SPACEFILL;
            else if (!strcmp(name, "licorice"))  r.type = MD_RENDER_MOL_REP_LICORICE;
            else if (!strcmp(name, "ball"))      r.type = MD_RENDER_MOL_REP_BALL_AND_STICK;
            else if (!strcmp(name, "bonds"))     r.type = MD_RENDER_MOL_REP_BONDS;
            else if (!strcmp(name, "cartoon"))   r.type = MD_RENDER_MOL_REP_CARTOON;
            else { fprintf(stderr, "unknown representation '%s'\n", name); return 2; }
            reps[n_reps++] = r;
        } else if ((!strcmp(a, "--fade") || !strcmp(a, "--fade-in")) && i + 3 < argc) {
            if (n_reps == 0) { fprintf(stderr, "%s before any --rep\n", a); return 2; }
            md_render_fade_t* f = &reps[n_reps - 1].fade;
            uint32_t mode = 0;
            if (!parse_fade_mode(argv[i + 1], &mode)) { fprintf(stderr, "fade mode: thin, shrink or both\n"); return 2; }
            f->mode |= mode;
            float* band = !strcmp(a, "--fade") ? f->out : f->in;
            band[0] = (float)atof(argv[i + 2]);
            band[1] = (float)atof(argv[i + 3]);
            i += 3;
        } else { fprintf(stderr, "unknown option '%s'\n", a); return 2; }
    }
    if (n_reps == 0) { memset(&reps[0], 0, sizeof(reps[0])); reps[0].type = MD_RENDER_MOL_REP_SPACEFILL; reps[0].param = -1.0f; n_reps = 1; }
    if (n_views == 0) views[n_views++] = "front";
    if (replicate < 1) replicate = 1;

    printf("mdlib %s\n", MD_BENCH_SOURCE_REV);

    /* The system */
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(16));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    md_str_t path = md_str_from_cstr(structure);
    const md_io_format_info_t* fmt = md_io_format_info(md_io_format_from_file(path));
    if (!fmt || !fmt->system_init || !fmt->system_init(&sys, &state, path, NULL)) {
        fprintf(stderr, "could not load '%s' as a system\n", structure);
        return 1;
    }
    if (coords) {
        md_str_t cpath = md_str_from_cstr(coords);
        const md_io_format_info_t* cf = md_io_format_info(md_io_format_from_file(cpath));
        md_system_state_t cs = { .alloc = alloc };
        if (!cf || !cf->state_init || !cf->state_init(&cs, cpath, sys.atom.count) || cs.num_atoms != sys.atom.count) {
            fprintf(stderr, "could not read coordinates for %zu atoms from '%s'\n", (size_t)sys.atom.count, coords);
            return 1;
        }
        state = cs;
    }
    if (!md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL)) { fprintf(stderr, "md_system_infer failed\n"); return 1; }
    const uint32_t n = (uint32_t)sys.atom.count;
    printf("%s: %u atoms, %zu bonds\n", structure, n, (size_t)(sys.bond.count + sys.bond.noncovalent_count));

    /* GPU */
    md_gpu_device_desc_t dd = {0};
    dd.label = "md_render_snap";
    md_gpu_device_t dev = md_gpu_device_create(&dd);
    if (!dev) { fprintf(stderr, "no GPU device: %s\n", md_gpu_last_error()); return 1; }
    md_gpu_device_info_t info;
    md_gpu_device_info(dev, &info);
    printf("device: %s (%s)\n", info.name, info.driver);
    md_gpu_stream_t s = md_gpu_stream_default(dev, MD_GPU_STREAM_GRAPHICS);
    md_render_t r = md_render_create(dev, NULL);
    md_render_mol_t mol = md_render_mol_create(dev, NULL);
    if (!s || !r || !mol) { fprintf(stderr, "md_render could not start\n"); return 1; }

    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    if (!md_render_mol_plan(&plan, &sys, NULL, temp)) { fprintf(stderr, "md_render_mol_plan failed\n"); return 1; }
    md_render_mol_topo_t topo;
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};
    if (!md_render_mol_write(s, &topo, &plan, md_gpu_malloc(s, MD_GPU_MEM_DEVICE, plan.size).gpu, content)) { fprintf(stderr, "md_render_mol_write failed\n"); return 1; }
    md_temp_end(temp);

    uint32_t* colors = malloc((size_t)n * 4);
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t c = 0;
        if (sys.atom.type_idx && sys.atom.type.color && (size_t)sys.atom.type_idx[i] < sys.atom.type.count) c = sys.atom.type.color[sys.atom.type_idx[i]];
        if (!c) c = md_atomic_number_cpk_color(md_atom_atomic_number(&sys.atom, i));
        colors[i] = c | 0xFF000000u;
    }
    uint32_t words = (n + 31) / 32;
    uint32_t* sel = calloc(words ? words : 1, 4);
    uint32_t n_sel = 0;
    for (int64_t i = sel_beg < 0 ? 0 : sel_beg; i < sel_end && i < (int64_t)n; ++i) { sel[i / 32] |= 1u << (i % 32); n_sel++; }

    md_gpu_addr_t P   = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, (size_t)n * 12).gpu;
    md_gpu_addr_t C   = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, (size_t)n * 4).gpu;
    md_gpu_addr_t S   = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, (size_t)(words ? words : 1) * 4).gpu;
    md_gpu_addr_t dyn = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    if (!P || !C || !S || !dyn) { fprintf(stderr, "out of GPU memory\n"); return 1; }
    md_gpu_upload(s, P, state.xyz, (size_t)n * 12);
    md_gpu_upload(s, C, colors, (size_t)n * 4);
    md_gpu_upload(s, S, sel, (size_t)(words ? words : 1) * 4);
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    if (n_sel) atoms.selection = S;
    if (!md_render_mol_update(mol, s, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, content)) { fprintf(stderr, "md_render_mol_update failed\n"); return 1; }

    /* Instances along the cell (or the bounds) */
    double lo[3] = { 1e30, 1e30, 1e30 }, hi[3] = { -1e30, -1e30, -1e30 };
    for (uint32_t i = 0; i < n; ++i) {
        const float* p = (const float*)&state.xyz[i];
        for (int k = 0; k < 3; ++k) { if (p[k] < lo[k]) lo[k] = p[k]; if (p[k] > hi[k]) hi[k] = p[k]; }
    }
    double A[3][3] = { { hi[0] - lo[0], 0, 0 }, { 0, hi[1] - lo[1], 0 }, { 0, 0, hi[2] - lo[2] } };
    if (md_unitcell_flags(&state.unitcell) & (MD_UNITCELL_ORTHO | MD_UNITCELL_TRICLINIC)) md_unitcell_A_extract_double(A, &state.unitcell);
    const uint32_t n_inst = replicate * replicate * replicate;
    float* xf = malloc((size_t)n_inst * 64);
    for (uint32_t ix = 0, k = 0; ix < replicate; ++ix) for (uint32_t iy = 0; iy < replicate; ++iy) for (uint32_t iz = 0; iz < replicate; ++iz, ++k) {
        float* m = xf + 16 * k;
        mat_ident(m);
        for (int d = 0; d < 3; ++d) m[12 + d] = (float)(ix * A[0][d] + iy * A[1][d] + iz * A[2][d]);
    }
    double ctr[3], rad = 0;
    for (int d = 0; d < 3; ++d) {
        double ext = (replicate - 1) * (A[0][d] + A[1][d] + A[2][d]);
        ctr[d] = 0.5 * (lo[d] + hi[d] + ext);
        double h = 0.5 * (hi[d] - lo[d] + fabs(ext));
        rad += h * h;
    }
    rad = sqrt(rad);
    printf("%u instance(s), %llu atoms drawn; bounding radius %.1f A\n", n_inst, (unsigned long long)n * n_inst, rad);

    /* The cartoon: the secondary structure (the state's, else computed), the
       profile per scale asked for, and its own colours: helix red, sheet
       yellow, coil grey (by the CA atoms) */
    md_gpu_addr_t SSB = 0, CC = 0;
    if (topo.backbone_count) {
        const md_secondary_structure_t* ss = md_system_state_secondary_structure(&state, &sys);
        if (!ss && md_system_state_backbone_compute(&state, &sys)) ss = md_system_state_secondary_structure(&state, &sys);
        const size_t nseg = sys.protein_backbone.segment.count;
        if (ss) {
            SSB = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, nseg * 4).gpu;
            md_gpu_upload(s, SSB, ss, nseg * 4);
        }
        uint32_t* cc = malloc((size_t)n * 4);
        for (uint32_t i = 0; i < n; ++i) cc[i] = colors[i];
        for (size_t g = 0; g < nseg; ++g) {
            int32_t ca = sys.protein_backbone.segment.atoms[g].ca;
            if (ca < 0 || (uint32_t)ca >= n) continue;
            int t = ss ? (int)ss[g] : 0;
            cc[ca] = (t >= MD_SECONDARY_STRUCTURE_HELIX_310 && t <= MD_SECONDARY_STRUCTURE_HELIX_PI) ? 0xFF4040E0u
                   : (t == MD_SECONDARY_STRUCTURE_BETA_SHEET || t == MD_SECONDARY_STRUCTURE_BETA_BRIDGE) ? 0xFF30C8F0u : 0xFFC8C8C8u;
        }
        CC = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, (size_t)n * 4).gpu;
        md_gpu_upload(s, CC, cc, (size_t)n * 4);
        free(cc);
        printf("cartoon: %u residues in %u pieces, secondary structure %s\n", topo.backbone_count, topo.backbone_chain_count, ss ? "given or computed" : "none (all coil)");
    }

    md_render_batch_t batches[64];
    uint32_t nb = 0;
    for (uint32_t k = 0; k < n_reps; ++k) {
        md_render_mol_rep_t rep;
        memset(&rep, 0, sizeof(rep));
        rep.key  = k + 1;
        rep.type = reps[k].type;
        float p = reps[k].param;
        switch (rep.type) {
        case MD_RENDER_MOL_REP_SPACEFILL:      rep.params.spacefill.radius_scale = p >= 0 ? p : 1.0f; break;
        case MD_RENDER_MOL_REP_LICORICE:       rep.params.licorice.radius = p >= 0 ? p : 0.2f; break;
        case MD_RENDER_MOL_REP_BALL_AND_STICK: rep.params.ball_and_stick.ball_scale = p >= 0 ? p : 0.3f; rep.params.ball_and_stick.stick_radius = 0.12f; break;
        case MD_RENDER_MOL_REP_BONDS:          rep.params.bonds.radius = p >= 0 ? p : 0.15f; break;
        case MD_RENDER_MOL_REP_CARTOON: {
            const float sc = p >= 0 ? p : 1.0f, scale[3] = { sc, sc, sc };
            md_gpu_addr_t prof = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, md_render_mol_cartoon_size(&topo) ? md_render_mol_cartoon_size(&topo) : 256).gpu;
            if (topo.backbone_count && !md_render_mol_cartoon_profile(mol, s, &topo, prof, SSB, scale, content)) { fprintf(stderr, "md_render_mol_cartoon_profile failed\n"); return 1; }
            rep.params.cartoon.profile = prof;
            break;
        }
        default: break;
        }
        rep.color.mode = MD_RENDER_COLOR_POINTS;
        rep.atom_color = rep.type == MD_RENDER_MOL_REP_CARTOON && CC ? CC : C;
        rep.fade = reps[k].fade;
        if (n_inst > 1) { rep.transform = xf; rep.transform_count = n_inst; }
        nb += md_render_mol_emit(&topo, dyn, &atoms, &rep, batches + nb, 64 - nb);
    }

    crc_init();
    uint8_t* rgb = malloc((size_t)W * H * 3);
    const float aspect = (float)W / (float)H;
    for (uint32_t vi = 0; vi < n_views; ++vi) {
        double dir[3];
        const char* vn = views[vi];
        if      (!strcmp(vn, "front")) { dir[0] = 0; dir[1] = 0; dir[2] = 1; }
        else if (!strcmp(vn, "back"))  { dir[0] = 0; dir[1] = 0; dir[2] = -1; }
        else if (!strcmp(vn, "side"))  { dir[0] = 1; dir[1] = 0; dir[2] = 0; }
        else if (!strcmp(vn, "top"))   { dir[0] = 0; dir[1] = 1; dir[2] = 0; }
        else if (!strcmp(vn, "iso"))   { dir[0] = 0.577; dir[1] = 0.577; dir[2] = 0.577; }
        else { fprintf(stderr, "unknown view '%s'\n", vn); continue; }
        const double fit = rad / sin(atan(0.4142 * (aspect < 1 ? aspect : 1.0)));
        const double dist = fit * zoom;
        md_render_frame_t f;
        memset(&f, 0, sizeof(f));
        f.width = W; f.height = H;
        f.camera = camera_look(ctr, dir, dist, aspect, dist > rad ? 0.05 * (dist - rad) + 0.01 : 0.05);
        f.batches = batches; f.batch_count = nb;
        f.jitter = MD_RENDER_JITTER_NONE; f.samples = 1;
        memcpy(f.content, content, sizeof(content));
        md_render_view_t v = md_render_view_create(r, "snap");
        md_tick_t t0 = md_tick_now();
        bool ok = md_render_view_draw(s, v, &f, NULL);
        md_gpu_stream_sync(s);
        double ms0 = md_tick_to_milliseconds(md_tick_now() - t0);
        /* again, against the first raster's pyramid (what a still view settles into) */
        f.content[MD_RENDER_CONTENT_GEOMETRY]++;
        t0 = md_tick_now();
        ok = ok && md_render_view_draw(s, v, &f, NULL);
        md_gpu_stream_sync(s);
        double ms1 = md_tick_to_milliseconds(md_tick_now() - t0);
        if (ok && timed) {
            /* `timed` draws after the geometry changed (the refit raster), best and
               median; then one with each pass synced and timed */
            double* ts = malloc(sizeof(double) * timed);
            for (uint32_t k = 0; k < timed; ++k) {
                f.content[MD_RENDER_CONTENT_GEOMETRY]++;
                md_tick_t tk = md_tick_now();
                ok = ok && md_render_view_draw(s, v, &f, NULL);
                md_gpu_stream_sync(s);
                ts[k] = md_tick_to_milliseconds(md_tick_now() - tk);
            }
            qsort(ts, timed, sizeof(double), cmp_double);
            printf("%s: %u draws after a geometry change: best %.2f, median %.2f ms\n", vn, timed, ts[0], ts[timed / 2]);
            free(ts);
            pass_prof_t prof = { .s = s };
            f.content[MD_RENDER_CONTENT_GEOMETRY]++;
            md_gpu_stream_sync(s);
            prof.t = md_tick_now();
            md_render_debug_pass_hook(pass_hook, &prof);
            ok = ok && md_render_view_draw(s, v, &f, NULL);
            md_render_debug_pass_hook(NULL, NULL);
            md_gpu_stream_sync(s);
            printf("  passes, each synced:");
            for (uint32_t k = 0; k < prof.n; ++k) printf(" %s %.2f", prof.name[k], prof.ms[k]);
            printf("\n");
        }
        md_render_stats_t st = {0};
        md_render_view_stats(v, &st);
        char file[1024];
        snprintf(file, sizeof(file), "%s_%s.png", out_prefix, vn);
        if (ok && compose(s, v, W, H, rgb) && write_png(file, rgb, W, H)) {
            printf("%s: draw %.2f ms, again %.2f ms; elements %llu: faded %u, outside %u, occluded %u, hardware %u, compute %u; scale at the centre %.2f px/A\n",
                   file, ms0, ms1, (unsigned long long)st.elements, st.faded, st.outside, st.occluded,
                   st.hardware[0] + st.hardware[1], st.compute[0] + st.compute[1], 0.5 * H / 0.4142 / dist);
        } else {
            fprintf(stderr, "%s: failed\n", file);
        }
        md_render_view_destroy(s, v);
    }
    md_gpu_stream_sync(s);
    free(rgb); free(xf); free(colors); free(sel);
    md_render_mol_destroy(mol);
    md_render_destroy(r);
    md_gpu_device_destroy(dev);
    md_system_free(&sys);
    md_vm_arena_destroy(alloc);
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: md_render_snap <structure> [--coords FILE] [--rep TYPE[:P]] [--fade MODE S0 S1] [--fade-in MODE S0 S1]\n"
                        "                      [--view front|back|side|top|iso] [--zoom Z] [--replicate N] [--select A:B]\n"
                        "                      [--res WxH] [--out PREFIX] [--time N]\n");
        return 2;
    }
    int rc = run(argc - 1, argv + 1);
    if (rc == 2) fprintf(stderr, "bad arguments (run without any for the usage)\n");
    return rc;
}

#else
int main(void) {
    fprintf(stderr, "md_render_snap needs MD_ENABLE_GPU\n");
    return 1;
}
#endif

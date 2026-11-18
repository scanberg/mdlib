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
//     --residues S0 S1     the last --rep (spacefill) hands its atoms over to a sphere per
//                          residue across S0 > S1 (pixels per A)
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

#include "render_scene.h"
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

/* ---- Run ----------------------------------------------------------------------- */

static int run(int argc, char** argv) {
    if (argc < 1) return 2;
    scene_desc_t sd;
    memset(&sd, 0, sizeof(sd));
    sd.structure = argv[0];
    sd.sel_beg = sd.sel_end = -1;
    const char* out_prefix = "snap";
    uint32_t W = 1280, H = 720, timed = 0;
    double zoom = 1.0;
    const char* views[8];
    uint32_t n_views = 0;
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (!strcmp(a, "--coords") && i + 1 < argc) sd.coords = argv[++i];
        else if (!strcmp(a, "--out") && i + 1 < argc) out_prefix = argv[++i];
        else if (!strcmp(a, "--res") && i + 1 < argc) { if (sscanf(argv[++i], "%ux%u", &W, &H) != 2) return 2; }
        else if (!strcmp(a, "--zoom") && i + 1 < argc) zoom = atof(argv[++i]);
        else if (!strcmp(a, "--replicate") && i + 1 < argc) sd.replicate = (uint32_t)atoi(argv[++i]);
        else if (!strcmp(a, "--time") && i + 1 < argc) timed = (uint32_t)atoi(argv[++i]);
        else if (!strcmp(a, "--select") && i + 1 < argc) { long long x, y; if (sscanf(argv[++i], "%lld:%lld", &x, &y) != 2) return 2; sd.sel_beg = x; sd.sel_end = y; }
        else if (!strcmp(a, "--gaps")) sd.chain_gaps = true;
        else if (!strcmp(a, "--orders")) sd.bond_orders = true;
        else if (!strcmp(a, "--view") && i + 1 < argc) { if (n_views < 8) views[n_views++] = argv[++i]; else ++i; }
        else if (!strcmp(a, "--rep") && i + 1 < argc) {
            if (sd.rep_count == 8) return 2;
            if (!scene_parse_rep(argv[++i], &sd.reps[sd.rep_count])) return 2;
            sd.rep_count++;
        } else if ((!strcmp(a, "--fade") || !strcmp(a, "--fade-in")) && i + 3 < argc) {
            if (sd.rep_count == 0) { fprintf(stderr, "%s before any --rep\n", a); return 2; }
            md_render_fade_t* f = &sd.reps[sd.rep_count - 1].fade;
            uint32_t mode = 0;
            if (!scene_parse_fade_mode(argv[i + 1], &mode)) { fprintf(stderr, "fade mode: thin, shrink or both\n"); return 2; }
            f->mode |= mode;
            float* band = !strcmp(a, "--fade") ? f->out : f->in;
            band[0] = (float)atof(argv[i + 2]);
            band[1] = (float)atof(argv[i + 3]);
            i += 3;
        } else if (!strcmp(a, "--residues") && i + 2 < argc) {
            if (sd.rep_count == 0) { fprintf(stderr, "%s before any --rep\n", a); return 2; }
            sd.reps[sd.rep_count - 1].residues[0] = (float)atof(argv[++i]);
            sd.reps[sd.rep_count - 1].residues[1] = (float)atof(argv[++i]);
        } else { fprintf(stderr, "unknown option '%s'\n", a); return 2; }
    }
    if (n_views == 0) views[n_views++] = "front";

    printf("mdlib %s\n", MD_BENCH_SOURCE_REV);

    scene_t sc;
    if (!scene_load(&sc, &sd)) return 1;

    /* GPU */
    md_gpu_device_desc_t dd = {0};
    dd.label = "md_render_snap";
    md_gpu_device_t dev = md_gpu_device_create(&dd);
    if (!dev) { fprintf(stderr, "no GPU device: %s\n", md_gpu_last_error()); return 1; }
    md_gpu_device_info_t info;
    md_gpu_device_info(dev, &info);
    printf("device: %s (%s)\n", info.name, info.driver);
    md_gpu_stream_t s = md_gpu_stream_default(dev, MD_GPU_STREAM_GRAPHICS);
    md_render_ctx_t r = md_render_create(dev, NULL);
    md_render_mol_ctx_t kernels = md_render_mol_create(dev, NULL);
    if (!s || !r || !kernels) { fprintf(stderr, "md_render could not start\n"); return 1; }
    if (!scene_upload(&sc, &sd, s, kernels)) return 1;
    const md_render_batch_t* batches = sc.batches;
    const uint32_t nb = sc.batch_count;

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
        const double dist = scene_fit_distance(&sc, aspect) * zoom;
        md_render_frame_t f;
        memset(&f, 0, sizeof(f));
        f.width = W; f.height = H;
        f.camera = camera_look(sc.ctr, dir, dist, aspect, scene_near(&sc, dist));
        f.batches = batches; f.batch_count = nb;
        f.jitter = MD_RENDER_JITTER_NONE; f.samples = 1;
        memcpy(f.content, sc.content, sizeof(sc.content));
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
    free(rgb);
    md_render_mol_destroy(kernels);
    md_render_destroy(r);
    md_gpu_device_destroy(dev);
    scene_free(&sc);
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: md_render_snap <structure> [--coords FILE] [--rep TYPE[:P]] [--fade MODE S0 S1] [--fade-in MODE S0 S1]\n"
                        "                      [--view front|back|side|top|iso] [--zoom Z] [--replicate N] [--select A:B]\n"
                        "                      [--res WxH] [--out PREFIX] [--time N] [--gaps] [--orders]\n");
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

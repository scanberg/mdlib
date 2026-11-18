// md_render_window: md_render_snap's scene in a window, presented through
// md_gpu's surface. For one question: does presentation work end to end on a
// machine, locally or remotely (ThinLinc, VNC, VirtualGL), and at what frame
// rate. It prints what it ran on and how it went; the picture is md_render_snap's
// look at the targets, not viamd's shading.
//
//   md_render_window <structure> [options]
//     --rep TYPE[:P]       as md_render_snap; repeatable (default spacefill)
//     --coords FILE        coordinates for it
//     --select A:B         atoms A .. B-1 selected (the mask)
//     --res WxH            window size (default 1280x720)
//     --still              the camera stands still: md_render draws once, then
//                          nothing is drawn or presented until the window
//                          changes (lazy rendering). Default: it turns
//     --turn DEG           degrees per second while turning (default 30)
//     --zoom Z             camera distance, 1 = the scene just fits (default 1)
//     --present MODE       vsync, mailbox or immediate (default vsync)
//     --format F           bgra8, srgb or rgba16f (default bgra8)
//     --adapter NAME       the GPU, as MD_GPU_DEVICE: an index or part of its name
//     --frames N           exit after N presented frames
//     --seconds S          exit after S seconds
//
// Keys: Escape closes, Space switches between turning and still.
//
// Every second: frames presented, their rate, the time between presents (mean
// and worst), and the host's time per frame split into recording and waiting
// (for a frame two back to finish, and for acquire: the presentation engine,
// which is what vsync paces). At the end: totals, and how many draws md_render
// found nothing to do for.

#include "render_scene.h"
#include <core/md_os.h>

#include "bench_rev.h"                    // MD_BENCH_SOURCE_REV

#if MD_ENABLE_GPU

#include <GLFW/glfw3.h>
#if defined(_WIN32)
#  define GLFW_EXPOSE_NATIVE_WIN32
#elif defined(__APPLE__)
#  define GLFW_EXPOSE_NATIVE_COCOA
#else
#  define GLFW_EXPOSE_NATIVE_X11
#endif
#include <GLFW/glfw3native.h>

#include "render_window_shaders.inl"

typedef struct compose_args_t {
    md_gpu_sampled_tex_t color, normal, mask;
    uint32_t             width, height;
    uint32_t             linear_out;
    uint32_t             _pad;
} compose_args_t;

static uint32_t g_fb_w, g_fb_h;
static bool     g_resized, g_refresh, g_toggle;

static void on_glfw_error(int code, const char* msg) { fprintf(stderr, "GLFW error 0x%x: %s\n", code, msg); }
static void on_fb_size(GLFWwindow* w, int width, int height) { (void)w; g_fb_w = (uint32_t)width; g_fb_h = (uint32_t)height; g_resized = true; }
static void on_refresh(GLFWwindow* w) { (void)w; g_refresh = true; }
static void on_key(GLFWwindow* w, int key, int scancode, int action, int mods) {
    (void)scancode; (void)mods;
    if (action != GLFW_PRESS) return;
    if (key == GLFW_KEY_ESCAPE) glfwSetWindowShouldClose(w, GLFW_TRUE);
    if (key == GLFW_KEY_SPACE)  g_toggle = true;
}

static const char* platform_name(int p) {
    switch (p) {
    case GLFW_PLATFORM_WIN32:   return "win32";
    case GLFW_PLATFORM_COCOA:   return "cocoa";
    case GLFW_PLATFORM_X11:     return "x11";
    case GLFW_PLATFORM_WAYLAND: return "wayland";
    case GLFW_PLATFORM_NULL:    return "null";
    default:                    return "unknown";
    }
}

static const char* format_name(md_gpu_format_t f) {
    switch (f) {
    case MD_GPU_FORMAT_BGRA8_UNORM:  return "BGRA8_UNORM";
    case MD_GPU_FORMAT_BGRA8_SRGB:   return "BGRA8_SRGB";
    case MD_GPU_FORMAT_RGBA16_FLOAT: return "RGBA16_FLOAT";
    default:                         return "other";
    }
}

static const char* present_name(md_gpu_present_mode_t m) {
    return m == MD_GPU_PRESENT_MAILBOX ? "mailbox" : m == MD_GPU_PRESENT_IMMEDIATE ? "immediate" : "vsync";
}

static int run(int argc, char** argv) {
    if (argc < 1) return 2;
    scene_desc_t sd;
    memset(&sd, 0, sizeof(sd));
    sd.structure = argv[0];
    sd.sel_beg = sd.sel_end = -1;
    uint32_t W = 1280, H = 720;
    uint64_t max_frames = 0;
    double max_seconds = 0, zoom = 1.0, turn = 30.0;
    bool still = false;
    const char* adapter = NULL;
    md_gpu_present_mode_t present = MD_GPU_PRESENT_VSYNC;
    md_gpu_format_t format = MD_GPU_FORMAT_BGRA8_UNORM;
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (!strcmp(a, "--coords") && i + 1 < argc) sd.coords = argv[++i];
        else if (!strcmp(a, "--res") && i + 1 < argc) { if (sscanf(argv[++i], "%ux%u", &W, &H) != 2) return 2; }
        else if (!strcmp(a, "--zoom") && i + 1 < argc) zoom = atof(argv[++i]);
        else if (!strcmp(a, "--turn") && i + 1 < argc) turn = atof(argv[++i]);
        else if (!strcmp(a, "--still")) still = true;
        else if (!strcmp(a, "--frames") && i + 1 < argc) max_frames = (uint64_t)atoll(argv[++i]);
        else if (!strcmp(a, "--seconds") && i + 1 < argc) max_seconds = atof(argv[++i]);
        else if (!strcmp(a, "--adapter") && i + 1 < argc) adapter = argv[++i];
        else if (!strcmp(a, "--select") && i + 1 < argc) { long long x, y; if (sscanf(argv[++i], "%lld:%lld", &x, &y) != 2) return 2; sd.sel_beg = x; sd.sel_end = y; }
        else if (!strcmp(a, "--present") && i + 1 < argc) {
            const char* m = argv[++i];
            if      (!strcmp(m, "vsync"))     present = MD_GPU_PRESENT_VSYNC;
            else if (!strcmp(m, "mailbox"))   present = MD_GPU_PRESENT_MAILBOX;
            else if (!strcmp(m, "immediate")) present = MD_GPU_PRESENT_IMMEDIATE;
            else return 2;
        } else if (!strcmp(a, "--format") && i + 1 < argc) {
            const char* f = argv[++i];
            if      (!strcmp(f, "bgra8"))   format = MD_GPU_FORMAT_BGRA8_UNORM;
            else if (!strcmp(f, "srgb"))    format = MD_GPU_FORMAT_BGRA8_SRGB;
            else if (!strcmp(f, "rgba16f")) format = MD_GPU_FORMAT_RGBA16_FLOAT;
            else return 2;
        } else if (!strcmp(a, "--rep") && i + 1 < argc) {
            if (sd.rep_count == 8) return 2;
            if (!scene_parse_rep(argv[++i], &sd.reps[sd.rep_count])) return 2;
            sd.rep_count++;
        } else { fprintf(stderr, "unknown option '%s'\n", a); return 2; }
    }

    printf("mdlib %s\n", MD_BENCH_SOURCE_REV);

    /* The window first: if there is no display, nothing else matters */
    glfwSetErrorCallback(on_glfw_error);
#if !defined(_WIN32) && !defined(__APPLE__)
    glfwInitHint(GLFW_PLATFORM, GLFW_PLATFORM_X11);         /* the one built here; Wayland desktops run it through XWayland */
#endif
    if (!glfwInit()) { fprintf(stderr, "glfwInit failed: no display?\n"); return 1; }
    printf("GLFW %s, platform %s\n", glfwGetVersionString(), platform_name(glfwGetPlatform()));
    glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
    GLFWwindow* win = glfwCreateWindow((int)W, (int)H, "md_render_window", NULL, NULL);
    if (!win) { fprintf(stderr, "could not create a window\n"); glfwTerminate(); return 1; }
    glfwSetFramebufferSizeCallback(win, on_fb_size);
    glfwSetWindowRefreshCallback(win, on_refresh);
    glfwSetKeyCallback(win, on_key);
    int fbw = 0, fbh = 0;
    glfwGetFramebufferSize(win, &fbw, &fbh);
    g_fb_w = (uint32_t)fbw; g_fb_h = (uint32_t)fbh;

    scene_t sc;
    if (!scene_load(&sc, &sd)) return 1;

    md_gpu_device_desc_t dd = {0};
    dd.label   = "md_render_window";
    dd.adapter = adapter;
    md_gpu_device_t dev = md_gpu_device_create(&dd);
    if (!dev) { fprintf(stderr, "no GPU device: %s\n", md_gpu_last_error()); return 1; }
    md_gpu_device_info_t info;
    md_gpu_device_info(dev, &info);
    printf("device: %s (%s), vendor 0x%04X, %s, present %s\n", info.name, info.driver, info.vendor_id,
           info.is_discrete ? "discrete" : "integrated or software", info.supports_present ? "supported" : "NOT supported");
    if (!info.supports_present) return 1;
    md_gpu_stream_t s = md_gpu_stream_default(dev, MD_GPU_STREAM_GRAPHICS);
    md_render_ctx_t r = md_render_create(dev, NULL);
    md_render_mol_ctx_t kernels = md_render_mol_create(dev, NULL);
    if (!s || !r || !kernels) { fprintf(stderr, "md_render could not start\n"); return 1; }
    if (!scene_upload(&sc, &sd, s, kernels)) return 1;

    /* The surface, from the window's native handles */
    md_gpu_surface_desc_t sfd = {0};
#if defined(_WIN32)
    sfd.system  = MD_GPU_WINDOW_WIN32;
    sfd.window  = (void*)glfwGetWin32Window(win);
#elif defined(__APPLE__)
    sfd.system  = MD_GPU_WINDOW_COCOA;
    sfd.window  = (void*)glfwGetCocoaView(win);
#else
    sfd.system  = MD_GPU_WINDOW_X11;
    sfd.window  = (void*)(uintptr_t)glfwGetX11Window(win);
    sfd.display = (void*)glfwGetX11Display();
#endif
    sfd.width        = g_fb_w;
    sfd.height       = g_fb_h;
    sfd.format       = format;
    sfd.present_mode = present;
    sfd.label        = "md_render_window";
    md_gpu_surface_t surf = md_gpu_surface_create(dev, &sfd);
    if (!surf && format != MD_GPU_FORMAT_BGRA8_UNORM) {
        fprintf(stderr, "surface in %s: %s; trying BGRA8_UNORM\n", format_name(format), md_gpu_last_error());
        sfd.format = MD_GPU_FORMAT_BGRA8_UNORM;
        surf = md_gpu_surface_create(dev, &sfd);
    }
    if (!surf) { fprintf(stderr, "no surface: %s\n", md_gpu_last_error()); return 1; }
    printf("surface: %ux%u, %s asked, %s asked (md_gpu falls back to vsync where it must)\n", g_fb_w, g_fb_h, format_name(sfd.format), present_name(present));

    md_render_view_t view = md_render_view_create(r, "window");
    md_gpu_pipeline_t pipe = NULL;
    md_gpu_format_t pipe_format = MD_GPU_FORMAT_INVALID;

    const double t_start = glfwGetTime();
    double t_last_present = 0, t_report = t_start, t_first = 0;
    uint64_t presented = 0, draws = 0, idle_draws = 0, acquire_fail = 0;
    uint32_t sec_frames = 0;
    double sec_gap_sum = 0, sec_gap_max = 0, sec_rec_sum = 0, sec_wait_sum = 0;
    bool printed_idle = false;
    md_gpu_sync_t frame_done[2] = { md_gpu_sync_none(), md_gpu_sync_none() };
    bool dirty = true;
    int rc = 0;

    while (!glfwWindowShouldClose(win)) {
        if (still && !dirty && !g_resized && !g_refresh) glfwWaitEventsTimeout(0.25);
        else glfwPollEvents();
        const double now = glfwGetTime();
        if (max_seconds > 0 && now - t_start >= max_seconds) break;
        if (g_toggle) { g_toggle = false; still = !still; printf("%s\n", still ? "still" : "turning"); dirty = true; }
        if (g_resized) { g_resized = false; md_gpu_surface_resize(surf, g_fb_w, g_fb_h); dirty = true; }
        if (g_refresh) { g_refresh = false; dirty = true; }
        if (g_fb_w == 0 || g_fb_h == 0) continue;           /* minimised */

        const double rec0 = glfwGetTime();
        const float aspect = (float)g_fb_w / (float)g_fb_h;
        const double ang = still ? 0.0 : (now - t_start) * turn * 3.14159265358979 / 180.0;
        double dir[3] = { sin(ang), 0.3, cos(ang) };
        const double dl = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
        for (int k = 0; k < 3; ++k) dir[k] /= dl;
        const double dist = scene_fit_distance(&sc, aspect) * zoom;
        md_render_frame_t f;
        memset(&f, 0, sizeof(f));
        f.width = g_fb_w; f.height = g_fb_h;
        f.camera = camera_look(sc.ctr, dir, dist, aspect, scene_near(&sc, dist));
        f.batches = sc.batches; f.batch_count = sc.batch_count;
        f.jitter = MD_RENDER_JITTER_NONE; f.samples = 1;
        f.time = now - t_start;
        memcpy(f.content, sc.content, sizeof(sc.content));

        double w0 = glfwGetTime();
        md_gpu_sync_wait(frame_done[presented % 2]);          /* at most two frames ahead */
        double waited = glfwGetTime() - w0;
        md_render_status_t st = {0};
        if (!md_render_view_draw(s, view, &f, &st)) { fprintf(stderr, "md_render_view_draw failed\n"); rc = 1; break; }
        draws++;
        if (st.work == 0 && !dirty) {                         /* nothing changed: nothing to present */
            idle_draws++;
            if (!printed_idle) { printf("idle after %llu presented frames: md_render had nothing to do, nothing presented\n", (unsigned long long)presented); printed_idle = true; }
            continue;
        }
        printed_idle = false;

        w0 = glfwGetTime();
        md_gpu_texture_t back = md_gpu_surface_acquire(s, surf);
        waited += glfwGetTime() - w0;
        if (!back) {
            if (++acquire_fail <= 3) fprintf(stderr, "acquire failed: %s\n", md_gpu_last_error());
            if (acquire_fail >= 10 && presented == 0) { fprintf(stderr, "giving up: no image could be acquired\n"); rc = 1; break; }
            continue;
        }
        const md_gpu_texture_desc_t* bd = md_gpu_texture_desc(back);
        if (!pipe || bd->format != pipe_format) {
            if (pipe) md_gpu_pipeline_destroy(pipe);
            md_gpu_pipeline_desc_t pd = {0};
            pd.vertex          = md_shader_render_window_full_vs_shader();
            pd.fragment        = md_shader_render_window_compose_fs_shader();
            pd.topology        = MD_GPU_TOPOLOGY_TRIANGLES;
            pd.color_count     = 1;
            pd.color[0].format = bd->format;
            pd.label           = "md_render_window.compose";
            pipe = md_gpu_pipeline_create(dev, &pd);
            pipe_format = bd->format;
            if (!pipe) { fprintf(stderr, "compose pipeline for %s: %s\n", format_name(bd->format), md_gpu_last_error()); rc = 1; break; }
            printf("swapchain image: %ux%u %s\n", bd->width, bd->height, format_name(bd->format));
        }
        md_render_targets_t t = md_render_view_targets(view);
        compose_args_t ca = {0};
        ca.color      = md_gpu_texture_sampled(t.color);
        ca.normal     = md_gpu_texture_sampled(t.normal);
        ca.mask       = md_gpu_texture_sampled(t.mask);
        ca.width      = g_fb_w;
        ca.height     = g_fb_h;
        ca.linear_out = bd->format == MD_GPU_FORMAT_BGRA8_SRGB || bd->format == MD_GPU_FORMAT_RGBA16_FLOAT;
        md_gpu_render_desc_t rd = {0};
        rd.color_count      = 1;
        rd.color[0].texture = back;
        rd.color[0].load    = MD_GPU_LOAD_DONT_CARE;
        bool ok = md_gpu_render_begin(s, &rd);
        ok = ok && md_gpu_draw(s, pipe, 3, 1, &ca, sizeof(ca));
        ok = md_gpu_render_end(s) && ok;
        ok = md_gpu_surface_present(s, surf) && ok;
        if (!ok) { fprintf(stderr, "compose or present failed: %s\n", md_gpu_last_error()); rc = 1; break; }
        frame_done[presented % 2] = md_gpu_stream_record(s);
        md_gpu_device_poll(dev);
        dirty = false;

        const double done = glfwGetTime();
        presented++;
        if (presented == 1) { t_first = done - t_start; printf("first frame presented %.1f ms after start\n", t_first * 1000.0); }
        else {
            const double gap = done - t_last_present;
            sec_gap_sum += gap;
            if (gap > sec_gap_max) sec_gap_max = gap;
        }
        sec_rec_sum  += done - rec0 - waited;
        sec_wait_sum += waited;
        sec_frames++;
        t_last_present = done;
        if (done - t_report >= 1.0) {
            const uint32_t gaps = sec_frames > 1 ? sec_frames - 1 : 1;
            printf("%5.1f s: %3u frames, %6.1f fps, between presents %6.2f ms mean %6.2f worst; host per frame: recording %5.2f ms, waiting %6.2f ms\n",
                   done - t_start, sec_frames, sec_frames / (done - t_report), 1000.0 * sec_gap_sum / gaps, 1000.0 * sec_gap_max,
                   1000.0 * sec_rec_sum / sec_frames, 1000.0 * sec_wait_sum / sec_frames);
            t_report = done; sec_frames = 0; sec_gap_sum = sec_gap_max = sec_rec_sum = sec_wait_sum = 0;
        }
        if (max_frames && presented >= max_frames) break;
    }

    md_gpu_stream_sync(s);
    const double total = glfwGetTime() - t_start;
    printf("done: %llu frames presented in %.1f s (%.1f fps overall), %llu draws, %llu with nothing to do, %llu failed acquires\n",
           (unsigned long long)presented, total, presented / (total > 0 ? total : 1), (unsigned long long)draws,
           (unsigned long long)idle_draws, (unsigned long long)acquire_fail);
    if (presented == 0) rc = 1;

    md_render_view_destroy(s, view);
    if (pipe) md_gpu_pipeline_destroy(pipe);
    md_gpu_surface_destroy(surf);
    md_gpu_stream_sync(s);
    md_render_mol_destroy(kernels);
    md_render_destroy(r);
    md_gpu_device_destroy(dev);
    scene_free(&sc);
    glfwDestroyWindow(win);
    glfwTerminate();
    return rc;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: md_render_window <structure> [--rep TYPE[:P]] [--coords FILE] [--select A:B] [--res WxH]\n"
                        "                        [--still] [--turn DEG] [--zoom Z] [--present vsync|mailbox|immediate]\n"
                        "                        [--format bgra8|srgb|rgba16f] [--adapter NAME] [--frames N] [--seconds S]\n");
        return 2;
    }
    int rc = run(argc - 1, argv + 1);
    if (rc == 2) fprintf(stderr, "bad arguments (run without any for the usage)\n");
    return rc;
}

#else
int main(void) {
    fprintf(stderr, "md_render_window needs MD_ENABLE_GPU\n");
    return 1;
}
#endif

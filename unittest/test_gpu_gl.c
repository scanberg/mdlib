#include "utest.h"

/* md_gpu_gl: textures md_gpu writes and an OpenGL context reads
   (MD_GPU_TEX_SHARED). Linux, where a GL context needs no window: EGL,
   surfaceless, loaded at run time, so the tests build without GL and skip
   where there is none. What runs where:

     memory      md_gpu's writes, read in GL, twice through the same memory
     semaphores  only where both APIs have them (not llvmpipe and lavapipe:
                 their texture is ordered on the host, as on Metal) */

#include <core/md_gpu.h>
#include <core/md_gpu_gl.h>

#include <stdint.h>
#include <string.h>

#if MD_ENABLE_GPU && defined(__linux__)

#include <dlfcn.h>
#include <unistd.h>

/* EGL, as much as a surfaceless core context needs */
typedef void* egl_display_t;
typedef void* egl_context_t;
typedef int   egl_int_t;
enum {
    EGL_NONE_                            = 0x3038,
    EGL_OPENGL_API_                      = 0x30A2,
    EGL_CONTEXT_MAJOR_VERSION_           = 0x3098,
    EGL_CONTEXT_MINOR_VERSION_           = 0x30FB,
    EGL_CONTEXT_OPENGL_PROFILE_MASK_     = 0x30FD,
    EGL_CONTEXT_OPENGL_CORE_PROFILE_BIT_ = 0x0001,
    EGL_PLATFORM_SURFACELESS_MESA_       = 0x31DD,
};

typedef struct egl_t {
    void*         lib;
    void*         (*GetProcAddress)(const char*);
    egl_display_t (*GetPlatformDisplayEXT)(unsigned, void*, const egl_int_t*);
    unsigned      (*Initialize)(egl_display_t, egl_int_t*, egl_int_t*);
    unsigned      (*BindAPI)(unsigned);
    egl_context_t (*CreateContext)(egl_display_t, void*, egl_context_t, const egl_int_t*);
    unsigned      (*MakeCurrent)(egl_display_t, void*, void*, egl_context_t);
    unsigned      (*DestroyContext)(egl_display_t, egl_context_t);
    unsigned      (*Terminate)(egl_display_t);
    egl_display_t display;
    egl_context_t context;
} egl_t;

static egl_t egl;

static void* gl_loader(const char* name) { return egl.GetProcAddress(name); }

/* A GL 4.5 core context current on this thread, or false with the reason */
static bool gl_open(const char** reason) {
    memset(&egl, 0, sizeof(egl));
    egl.lib = dlopen("libEGL.so.1", RTLD_NOW | RTLD_LOCAL);
    if (!egl.lib) { *reason = "no libEGL.so.1"; return false; }
    *(void**)&egl.GetProcAddress = dlsym(egl.lib, "eglGetProcAddress");
    if (!egl.GetProcAddress) { *reason = "no eglGetProcAddress"; return false; }
    *(void**)&egl.GetPlatformDisplayEXT = egl.GetProcAddress("eglGetPlatformDisplayEXT");
    *(void**)&egl.Initialize     = dlsym(egl.lib, "eglInitialize");
    *(void**)&egl.BindAPI        = dlsym(egl.lib, "eglBindAPI");
    *(void**)&egl.CreateContext  = dlsym(egl.lib, "eglCreateContext");
    *(void**)&egl.MakeCurrent    = dlsym(egl.lib, "eglMakeCurrent");
    *(void**)&egl.DestroyContext = dlsym(egl.lib, "eglDestroyContext");
    *(void**)&egl.Terminate      = dlsym(egl.lib, "eglTerminate");
    if (!egl.GetPlatformDisplayEXT || !egl.Initialize || !egl.CreateContext) { *reason = "EGL lacks surfaceless displays"; return false; }
    egl.display = egl.GetPlatformDisplayEXT(EGL_PLATFORM_SURFACELESS_MESA_, NULL, NULL);
    if (!egl.display || !egl.Initialize(egl.display, NULL, NULL)) { *reason = "no surfaceless EGL display"; return false; }
    egl.BindAPI(EGL_OPENGL_API_);
    const egl_int_t attr[] = { EGL_CONTEXT_MAJOR_VERSION_, 4, EGL_CONTEXT_MINOR_VERSION_, 5,
                               EGL_CONTEXT_OPENGL_PROFILE_MASK_, EGL_CONTEXT_OPENGL_CORE_PROFILE_BIT_, EGL_NONE_ };
    egl.context = egl.CreateContext(egl.display, NULL, NULL, attr);
    if (!egl.context || !egl.MakeCurrent(egl.display, NULL, NULL, egl.context)) { *reason = "no GL 4.5 core context"; return false; }
    return true;
}

static void gl_close(void) {
    if (egl.context) { egl.MakeCurrent(egl.display, NULL, NULL, NULL); egl.DestroyContext(egl.display, egl.context); }
    if (egl.display) egl.Terminate(egl.display);
    if (egl.lib) dlclose(egl.lib);
    memset(&egl, 0, sizeof(egl));
}

/* The texels of a GL texture, through glGetTexImage */
static bool gl_read(const md_gpu_gl_texture_t* t, unsigned format, unsigned type, void* out) {
    void (*BindTexture)(unsigned, unsigned) = (void (*)(unsigned, unsigned))gl_loader("glBindTexture");
    void (*GetTexImage)(unsigned, int, unsigned, unsigned, void*) = (void (*)(unsigned, int, unsigned, unsigned, void*))gl_loader("glGetTexImage");
    unsigned (*GetError)(void) = (unsigned (*)(void))gl_loader("glGetError");
    if (!BindTexture || !GetTexImage || !GetError) return false;
    BindTexture(t->target, t->name);
    GetTexImage(t->target, 0, format, type, out);
    BindTexture(t->target, 0);
    return GetError() == 0;
}

typedef struct gpu_gl_t {
    md_gpu_device_t dev;
    md_gpu_stream_t s;
} gpu_gl_t;

/* A device that shares, and a GL context; or a reason to skip */
static const char* gpu_gl_open(gpu_gl_t* f) {
    memset(f, 0, sizeof(*f));
    const char* reason = NULL;
    if (!gl_open(&reason)) { gl_close(); return reason; }
    f->dev = md_gpu_device_create(NULL);
    if (!f->dev) { gl_close(); return "no GPU device"; }
    md_gpu_device_info_t info;
    md_gpu_device_info(f->dev, &info);
    if (!info.supports_sharing) { md_gpu_device_destroy(f->dev); gl_close(); return "the device cannot share textures"; }
    f->s = md_gpu_stream_default(f->dev, MD_GPU_STREAM_COMPUTE);
    return NULL;
}

static void gpu_gl_close(gpu_gl_t* f) {
    if (f->dev) md_gpu_device_destroy(f->dev);
    gl_close();
}

/* What md_gpu writes, GL reads: the same memory, round after round */
UTEST(gpu_gl, shared_texture_reads_in_gl) {
    gpu_gl_t f;
    const char* skip = gpu_gl_open(&f);
    if (skip) UTEST_SKIP(skip);

    enum { W = 37, H = 21 };
    md_gpu_texture_desc_t td = {0};
    td.type   = MD_GPU_TEX_2D;
    td.format = MD_GPU_FORMAT_RGBA8_UNORM;
    td.usage  = MD_GPU_TEX_SAMPLED | MD_GPU_TEX_SHARED;
    td.width  = W;
    td.height = H;
    td.label  = "shared";
    md_gpu_texture_t tex = md_gpu_texture_create(f.s, &td);
    ASSERT_TRUE(tex != NULL);

    md_gpu_gl_texture_t gt;
    ASSERT_TRUE(md_gpu_gl_import(&gt, tex, gl_loader));
    EXPECT_EQ(0x0DE1u, gt.target);   /* GL_TEXTURE_2D */

    static uint32_t src[W * H], got[W * H];
    for (int round = 0; round < 3; ++round) {
        for (uint32_t i = 0; i < W * H; ++i) src[i] = (i * 2654435761u) ^ (uint32_t)round * 0x01010101u;
        ASSERT_TRUE(md_gpu_texture_share_acquire(f.s, tex));
        ASSERT_TRUE(md_gpu_upload_texture(f.s, tex, NULL, src, sizeof(src)));
        ASSERT_TRUE(md_gpu_texture_share_release(f.s, tex));

        md_gpu_gl_begin_read(&gt);
        memset(got, 0, sizeof(got));
        ASSERT_TRUE(gl_read(&gt, 0x1908 /* GL_RGBA */, 0x1401 /* GL_UNSIGNED_BYTE */, got));
        md_gpu_gl_end_read(&gt);
        uint32_t wrong = 0;
        for (uint32_t i = 0; i < W * H; ++i) wrong += got[i] != src[i];
        EXPECT_EQ(0u, wrong);
    }

    md_gpu_gl_release(&gt);
    md_gpu_texture_destroy(tex);
    md_gpu_stream_sync(f.s);
    gpu_gl_close(&f);
}

/* A float target too (depth goes over as R32_FLOAT) */
UTEST(gpu_gl, shared_float_texture) {
    gpu_gl_t f;
    const char* skip = gpu_gl_open(&f);
    if (skip) UTEST_SKIP(skip);

    enum { W = 16, H = 9 };
    md_gpu_texture_desc_t td = {0};
    td.type   = MD_GPU_TEX_2D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.usage  = MD_GPU_TEX_STORAGE | MD_GPU_TEX_SHARED;
    td.width  = W;
    td.height = H;
    md_gpu_texture_t tex = md_gpu_texture_create(f.s, &td);
    ASSERT_TRUE(tex != NULL);
    md_gpu_gl_texture_t gt;
    ASSERT_TRUE(md_gpu_gl_import(&gt, tex, gl_loader));

    static float src[W * H], got[W * H];
    for (uint32_t i = 0; i < W * H; ++i) src[i] = 0.25f + (float)i / 7.0f;
    ASSERT_TRUE(md_gpu_texture_share_acquire(f.s, tex));
    ASSERT_TRUE(md_gpu_upload_texture(f.s, tex, NULL, src, sizeof(src)));
    ASSERT_TRUE(md_gpu_texture_share_release(f.s, tex));
    md_gpu_gl_begin_read(&gt);
    ASSERT_TRUE(gl_read(&gt, 0x1903 /* GL_RED */, 0x1406 /* GL_FLOAT */, got));
    md_gpu_gl_end_read(&gt);
    EXPECT_EQ(0, memcmp(src, got, sizeof(src)));

    md_gpu_gl_release(&gt);
    md_gpu_texture_destroy(tex);
    md_gpu_stream_sync(f.s);
    gpu_gl_close(&f);
}

/* The hand-over's rules, by name: no second release before an acquire, no
   sharing a texture made without SHARED, nor placing a shared one */
UTEST(gpu_gl, share_rules) {
    md_gpu_device_t dev = md_gpu_device_create(NULL);
    if (!dev) UTEST_SKIP("no GPU device");
    md_gpu_device_info_t info;
    md_gpu_device_info(dev, &info);
    if (!info.supports_sharing) { md_gpu_device_destroy(dev); UTEST_SKIP("the device cannot share textures"); }
    md_gpu_stream_t s = md_gpu_stream_default(dev, MD_GPU_STREAM_COMPUTE);

    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_2D; td.format = MD_GPU_FORMAT_RGBA8_UNORM; td.width = 8; td.height = 8;
    td.usage = MD_GPU_TEX_SAMPLED;
    md_gpu_texture_t plain = md_gpu_texture_create(s, &td);
    ASSERT_TRUE(plain != NULL);
    md_gpu_texture_share_t sh;
    EXPECT_FALSE(md_gpu_texture_share(plain, &sh));
    EXPECT_FALSE(md_gpu_texture_share_release(s, plain));

    td.usage = MD_GPU_TEX_SAMPLED | MD_GPU_TEX_SHARED;
    td.mip_levels = 2;
    EXPECT_TRUE(md_gpu_texture_create(s, &td) == NULL);           /* one mip */
    td.mip_levels = 1;
    td.format = MD_GPU_FORMAT_RG8_UNORM;
    EXPECT_TRUE(md_gpu_texture_create(s, &td) == NULL);           /* not a shared format */
    td.format = MD_GPU_FORMAT_RGBA8_UNORM;
    md_gpu_texture_t shared = md_gpu_texture_create(s, &td);
    ASSERT_TRUE(shared != NULL);
    ASSERT_TRUE(md_gpu_texture_share(shared, &sh));
    EXPECT_EQ(8u, sh.width);
    EXPECT_TRUE(sh.memory_size >= 8u * 8u * 4u);
    EXPECT_EQ(0, memcmp(sh.device_uuid, info.device_uuid, 16));
    if (sh.kind == MD_GPU_SHARE_FD) {                              /* ours: no import takes them over */
        close((int)sh.memory);
        if (sh.ready) { close((int)sh.ready); close((int)sh.done); }
    }

    if (!sh.ready) {   /* without semaphores: release waits on the host, and a second one is an error */
        EXPECT_TRUE(md_gpu_texture_share_release(s, shared));
        EXPECT_FALSE(md_gpu_texture_share_release(s, shared));
        EXPECT_TRUE(md_gpu_texture_share_acquire(s, shared));
    }

    md_gpu_texture_destroy(plain);
    md_gpu_texture_destroy(shared);
    md_gpu_stream_sync(s);
    md_gpu_device_destroy(dev);
}

#endif

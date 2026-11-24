// md_gpu_gl: see md_gpu_gl.h. The GL types and enums it needs, spelled out here, so that it needs no GL headers.

#include "md_gpu_gl.h"

#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_str.h>

typedef unsigned int   gl_enum_t;
typedef unsigned int   gl_uint_t;
typedef int            gl_int_t;
typedef int            gl_sizei_t;
typedef unsigned char  gl_ubyte_t;
typedef uint64_t       gl_uint64_t;

#if defined(_WIN32) && !defined(_WIN64)
#  define MD_GL_API __stdcall
#else
#  define MD_GL_API
#endif

enum {
    GL_NO_ERROR_                    = 0,
    GL_TEXTURE_2D_                  = 0x0DE1,
    GL_TEXTURE_RECTANGLE_           = 0x84F5,
    GL_TEXTURE_BINDING_2D_          = 0x8069,
    GL_TEXTURE_BINDING_RECTANGLE_   = 0x84F6,
    GL_TEXTURE_MAG_FILTER_          = 0x2800,
    GL_TEXTURE_MIN_FILTER_          = 0x2801,
    GL_TEXTURE_WRAP_S_              = 0x2802,
    GL_TEXTURE_WRAP_T_              = 0x2803,
    GL_NEAREST_                     = 0x2600,
    GL_CLAMP_TO_EDGE_               = 0x812F,
    GL_NUM_EXTENSIONS_              = 0x821D,
    GL_EXTENSIONS_                  = 0x1F03,
    GL_RED_                         = 0x1903,
    GL_RG_                          = 0x8227,
    GL_RGBA_                        = 0x1908,
    GL_UNSIGNED_BYTE_               = 0x1401,
    GL_HALF_FLOAT_                  = 0x140B,
    GL_FLOAT_                       = 0x1406,
    GL_RGBA8_                       = 0x8058,
    GL_SRGB8_ALPHA8_                = 0x8C43,
    GL_RGBA16F_                     = 0x881A,
    GL_RG16F_                       = 0x822F,
    GL_R32F_                        = 0x822E,
    GL_TEXTURE_TILING_EXT_          = 0x9580,
    GL_DEDICATED_MEMORY_OBJECT_EXT_ = 0x9581,
    GL_OPTIMAL_TILING_EXT_          = 0x9584,
    GL_HANDLE_TYPE_OPAQUE_FD_EXT_   = 0x9586,
    GL_HANDLE_TYPE_OPAQUE_WIN32_EXT_= 0x9587,
    GL_LAYOUT_GENERAL_EXT_          = 0x958D,
    GL_NUM_DEVICE_UUIDS_EXT_        = 0x9596,
    GL_DEVICE_UUID_EXT_             = 0x9597,
    GL_DRIVER_UUID_EXT_             = 0x9598,
};

typedef struct md_gl_fn_t {
    gl_enum_t        (MD_GL_API *GetError)(void);
    void             (MD_GL_API *GetIntegerv)(gl_enum_t, gl_int_t*);
    const gl_ubyte_t*(MD_GL_API *GetStringi)(gl_enum_t, gl_uint_t);
    void             (MD_GL_API *GenTextures)(gl_sizei_t, gl_uint_t*);
    void             (MD_GL_API *DeleteTextures)(gl_sizei_t, const gl_uint_t*);
    void             (MD_GL_API *BindTexture)(gl_enum_t, gl_uint_t);
    void             (MD_GL_API *TexParameteri)(gl_enum_t, gl_enum_t, gl_int_t);
    void             (MD_GL_API *Flush)(void);
    void             (MD_GL_API *Finish)(void);
    /* GL_EXT_memory_object, _fd, _win32 */
    void             (MD_GL_API *GetUnsignedBytevEXT)(gl_enum_t, gl_ubyte_t*);
    void             (MD_GL_API *GetUnsignedBytei_vEXT)(gl_enum_t, gl_uint_t, gl_ubyte_t*);
    void             (MD_GL_API *CreateMemoryObjectsEXT)(gl_sizei_t, gl_uint_t*);
    void             (MD_GL_API *DeleteMemoryObjectsEXT)(gl_sizei_t, const gl_uint_t*);
    void             (MD_GL_API *MemoryObjectParameterivEXT)(gl_uint_t, gl_enum_t, const gl_int_t*);
    void             (MD_GL_API *ImportMemoryFdEXT)(gl_uint_t, gl_uint64_t, gl_enum_t, gl_int_t);
    void             (MD_GL_API *ImportMemoryWin32HandleEXT)(gl_uint_t, gl_uint64_t, gl_enum_t, void*);
    void             (MD_GL_API *TexStorageMem2DEXT)(gl_enum_t, gl_sizei_t, gl_enum_t, gl_sizei_t, gl_sizei_t, gl_uint_t, gl_uint64_t);
    /* GL_EXT_semaphore, _fd, _win32 */
    void             (MD_GL_API *GenSemaphoresEXT)(gl_sizei_t, gl_uint_t*);
    void             (MD_GL_API *DeleteSemaphoresEXT)(gl_sizei_t, const gl_uint_t*);
    void             (MD_GL_API *ImportSemaphoreFdEXT)(gl_uint_t, gl_enum_t, gl_int_t);
    void             (MD_GL_API *ImportSemaphoreWin32HandleEXT)(gl_uint_t, gl_enum_t, void*);
    void             (MD_GL_API *WaitSemaphoreEXT)(gl_uint_t, gl_uint_t, const gl_uint_t*, gl_uint_t, const gl_uint_t*, const gl_enum_t*);
    void             (MD_GL_API *SignalSemaphoreEXT)(gl_uint_t, gl_uint_t, const gl_uint_t*, gl_uint_t, const gl_uint_t*, const gl_enum_t*);
    /* CGL (macOS) */
    void*            (*CGLGetCurrentContext)(void);
    int              (*CGLTexImageIOSurface2D)(void*, gl_enum_t, gl_enum_t, gl_sizei_t, gl_sizei_t, gl_enum_t, gl_enum_t, void*, gl_uint_t);
} md_gl_fn_t;

/* One table for the process: entry points do not depend on the context on
   the platforms this serves */
static md_gl_fn_t md_gl_fn;

#define MD_GL_LOAD(field, name) (*(void**)&md_gl_fn.field = load(name))

static void md_gl_load(md_gpu_gl_loader_t load) {
    MD_GL_LOAD(GetError,                      "glGetError");
    MD_GL_LOAD(GetIntegerv,                   "glGetIntegerv");
    MD_GL_LOAD(GetStringi,                    "glGetStringi");
    MD_GL_LOAD(GenTextures,                   "glGenTextures");
    MD_GL_LOAD(DeleteTextures,                "glDeleteTextures");
    MD_GL_LOAD(BindTexture,                   "glBindTexture");
    MD_GL_LOAD(TexParameteri,                 "glTexParameteri");
    MD_GL_LOAD(Flush,                         "glFlush");
    MD_GL_LOAD(Finish,                        "glFinish");
    MD_GL_LOAD(GetUnsignedBytevEXT,           "glGetUnsignedBytevEXT");
    MD_GL_LOAD(GetUnsignedBytei_vEXT,         "glGetUnsignedBytei_vEXT");
    MD_GL_LOAD(CreateMemoryObjectsEXT,        "glCreateMemoryObjectsEXT");
    MD_GL_LOAD(DeleteMemoryObjectsEXT,        "glDeleteMemoryObjectsEXT");
    MD_GL_LOAD(MemoryObjectParameterivEXT,    "glMemoryObjectParameterivEXT");
    MD_GL_LOAD(ImportMemoryFdEXT,             "glImportMemoryFdEXT");
    MD_GL_LOAD(ImportMemoryWin32HandleEXT,    "glImportMemoryWin32HandleEXT");
    MD_GL_LOAD(TexStorageMem2DEXT,            "glTexStorageMem2DEXT");
    MD_GL_LOAD(GenSemaphoresEXT,              "glGenSemaphoresEXT");
    MD_GL_LOAD(DeleteSemaphoresEXT,           "glDeleteSemaphoresEXT");
    MD_GL_LOAD(ImportSemaphoreFdEXT,          "glImportSemaphoreFdEXT");
    MD_GL_LOAD(ImportSemaphoreWin32HandleEXT, "glImportSemaphoreWin32HandleEXT");
    MD_GL_LOAD(WaitSemaphoreEXT,              "glWaitSemaphoreEXT");
    MD_GL_LOAD(SignalSemaphoreEXT,            "glSignalSemaphoreEXT");
    MD_GL_LOAD(CGLGetCurrentContext,          "CGLGetCurrentContext");
    MD_GL_LOAD(CGLTexImageIOSurface2D,        "CGLTexImageIOSurface2D");
}

/* Whether the context lists the extension. A loader can hand out an entry
   point for any name (Mesa's does), so a pointer proves nothing. */
static bool md_gl_has(const char* name) {
    if (!md_gl_fn.GetIntegerv || !md_gl_fn.GetStringi) return false;
    gl_int_t n = 0;
    md_gl_fn.GetIntegerv(GL_NUM_EXTENSIONS_, &n);
    const md_str_t want = md_str_from_cstr(name);
    for (gl_int_t i = 0; i < n; ++i) {
        const char* e = (const char*)md_gl_fn.GetStringi(GL_EXTENSIONS_, (gl_uint_t)i);
        if (e && md_str_eq(md_str_from_cstr(e), want)) return true;
    }
    return false;
}

/* The GL internal format, and for an IOSurface the format and type of its texels */
static bool md_gl_format(md_gpu_format_t f, gl_enum_t* internal, gl_enum_t* format, gl_enum_t* type) {
    switch (f) {
    case MD_GPU_FORMAT_RGBA8_UNORM:  *internal = GL_RGBA8_;        *format = GL_RGBA_; *type = GL_UNSIGNED_BYTE_; return true;
    case MD_GPU_FORMAT_RGBA8_SRGB:   *internal = GL_SRGB8_ALPHA8_; *format = GL_RGBA_; *type = GL_UNSIGNED_BYTE_; return true;
    case MD_GPU_FORMAT_RGBA16_FLOAT: *internal = GL_RGBA16F_;      *format = GL_RGBA_; *type = GL_HALF_FLOAT_;    return true;
    case MD_GPU_FORMAT_RG16_FLOAT:   *internal = GL_RG16F_;        *format = GL_RG_;   *type = GL_HALF_FLOAT_;    return true;
    case MD_GPU_FORMAT_R32_FLOAT:    *internal = GL_R32F_;         *format = GL_RED_;  *type = GL_FLOAT_;         return true;
    default: return false;
    }
}

static bool md_gl_same_gpu(const md_gpu_texture_share_t* sh) {
    if (!md_gl_fn.GetUnsignedBytevEXT || !md_gl_fn.GetUnsignedBytei_vEXT) return false;
    gl_ubyte_t driver[16] = {0};
    md_gl_fn.GetUnsignedBytevEXT(GL_DRIVER_UUID_EXT_, driver);
    if (md_memcmp(driver, sh->driver_uuid, 16) != 0) return false;
    gl_int_t n = 0;
    md_gl_fn.GetIntegerv(GL_NUM_DEVICE_UUIDS_EXT_, &n);
    for (gl_int_t i = 0; i < n; ++i) {
        gl_ubyte_t device[16] = {0};
        md_gl_fn.GetUnsignedBytei_vEXT(GL_DEVICE_UUID_EXT_, (gl_uint_t)i, device);
        if (md_memcmp(device, sh->device_uuid, 16) == 0) return true;
    }
    return false;
}

/* The texture parameters a read of it needs, with the 2D binding put back */
static void md_gl_texture_params(gl_enum_t target) {
    md_gl_fn.TexParameteri(target, GL_TEXTURE_MIN_FILTER_, GL_NEAREST_);
    md_gl_fn.TexParameteri(target, GL_TEXTURE_MAG_FILTER_, GL_NEAREST_);
    md_gl_fn.TexParameteri(target, GL_TEXTURE_WRAP_S_, GL_CLAMP_TO_EDGE_);
    md_gl_fn.TexParameteri(target, GL_TEXTURE_WRAP_T_, GL_CLAMP_TO_EDGE_);
}

static bool md_gl_import_iosurface(md_gpu_gl_texture_t* out, const md_gpu_texture_share_t* sh, gl_enum_t internal, gl_enum_t format, gl_enum_t type) {
    if (!md_gl_fn.CGLGetCurrentContext || !md_gl_fn.CGLTexImageIOSurface2D) {
        MD_LOG_ERROR("md_gpu_gl_import: the loader has no CGLTexImageIOSurface2D (OpenGL.framework)");
        return false;
    }
    void* ctx = md_gl_fn.CGLGetCurrentContext();
    if (!ctx) { MD_LOG_ERROR("md_gpu_gl_import: no current CGL context"); return false; }
    gl_int_t prev = 0;
    md_gl_fn.GetIntegerv(GL_TEXTURE_BINDING_RECTANGLE_, &prev);
    md_gl_fn.GenTextures(1, &out->name);
    md_gl_fn.BindTexture(GL_TEXTURE_RECTANGLE_, out->name);
    const int err = md_gl_fn.CGLTexImageIOSurface2D(ctx, GL_TEXTURE_RECTANGLE_, internal, (gl_sizei_t)sh->width, (gl_sizei_t)sh->height,
                                                    format, type, (void*)(uintptr_t)sh->memory, 0);
    md_gl_texture_params(GL_TEXTURE_RECTANGLE_);
    md_gl_fn.BindTexture(GL_TEXTURE_RECTANGLE_, (gl_uint_t)prev);
    if (err != 0) {
        MD_LOG_ERROR("md_gpu_gl_import: CGLTexImageIOSurface2D failed (%d)", err);
        md_gl_fn.DeleteTextures(1, &out->name);
        out->name = 0;
        return false;
    }
    out->target = GL_TEXTURE_RECTANGLE_;
    return true;
}

bool md_gpu_gl_import(md_gpu_gl_texture_t* out, md_gpu_texture_t tex, md_gpu_gl_loader_t load) {
    if (!out || !tex || !load) { MD_LOG_ERROR("md_gpu_gl_import: null argument"); return false; }
    md_memset(out, 0, sizeof(*out));
    md_gl_load(load);
    out->fn = &md_gl_fn;
    if (!md_gl_fn.GetError || !md_gl_fn.GetIntegerv || !md_gl_fn.GenTextures || !md_gl_fn.BindTexture || !md_gl_fn.TexParameteri ||
        !md_gl_fn.Flush || !md_gl_fn.Finish || !md_gl_fn.DeleteTextures) {
        MD_LOG_ERROR("md_gpu_gl_import: the loader lacks core GL entry points (no current context?)");
        return false;
    }

    const md_gpu_texture_desc_t* desc = md_gpu_texture_desc(tex);
    gl_enum_t internal = 0, format = 0, type = 0;
    if (!desc || !(desc->usage & MD_GPU_TEX_SHARED) || !md_gl_format(desc->format, &internal, &format, &type)) {
        MD_LOG_ERROR("md_gpu_gl_import: not a texture md_gpu shares (MD_GPU_TEX_SHARED, one of its formats)");
        return false;
    }

    /* What GL needs for the platform's kind, checked before md_gpu makes the
       handles to hand over (fds, which the import would take over) */
#if defined(__APPLE__)
    const md_gpu_share_kind_t kind = MD_GPU_SHARE_IOSURFACE;
#elif defined(_WIN32)
    const md_gpu_share_kind_t kind = MD_GPU_SHARE_WIN32;
#else
    const md_gpu_share_kind_t kind = MD_GPU_SHARE_FD;
#endif
    if (kind != MD_GPU_SHARE_IOSURFACE) {
        const char* mem_ext = kind == MD_GPU_SHARE_WIN32 ? "GL_EXT_memory_object_win32" : "GL_EXT_memory_object_fd";
        if (!md_gl_has("GL_EXT_memory_object") || !md_gl_has(mem_ext)) {
            MD_LOG_ERROR("md_gpu_gl_import: this GL context lacks GL_EXT_memory_object or %s", mem_ext);
            return false;
        }
    }

    md_gpu_texture_share_t sh;
    if (!md_gpu_texture_share(tex, &sh)) return false;
    out->kind = sh.kind;
    if (sh.kind == MD_GPU_SHARE_IOSURFACE) return md_gl_import_iosurface(out, &sh, internal, format, type);

    if (!md_gl_same_gpu(&sh)) {
        MD_LOG_ERROR("md_gpu_gl_import: the GL context runs on another GPU, or driver, than the md_gpu device "
                     "(GL_DEVICE_UUID_EXT / GL_DRIVER_UUID_EXT differ from md_gpu_device_info_t's)");
        return false;
    }
    const bool win32 = sh.kind == MD_GPU_SHARE_WIN32;
    if (sh.ready) {
        const char* sem_ext = win32 ? "GL_EXT_semaphore_win32" : "GL_EXT_semaphore_fd";
        if (!md_gl_has("GL_EXT_semaphore") || !md_gl_has(sem_ext)) {
            /* md_gpu would signal a `ready` GL never waits on, and wait for a `done` GL never signals */
            MD_LOG_ERROR("md_gpu_gl_import: the md_gpu device orders shared textures by semaphores, but this GL context lacks "
                         "GL_EXT_semaphore or %s", sem_ext);
            return false;
        }
    }
    while (md_gl_fn.GetError() != GL_NO_ERROR_) {}   /* errors before ours are not ours */

    /* The memory, dedicated as md_gpu allocated it, and the texture over it in
       the tiling md_gpu laid it out in */
    md_gl_fn.CreateMemoryObjectsEXT(1, &out->memory);
    const gl_int_t dedicated = 1;
    md_gl_fn.MemoryObjectParameterivEXT(out->memory, GL_DEDICATED_MEMORY_OBJECT_EXT_, &dedicated);
    if (win32) md_gl_fn.ImportMemoryWin32HandleEXT(out->memory, sh.memory_size, GL_HANDLE_TYPE_OPAQUE_WIN32_EXT_, (void*)(uintptr_t)sh.memory);
    else       md_gl_fn.ImportMemoryFdEXT(out->memory, sh.memory_size, GL_HANDLE_TYPE_OPAQUE_FD_EXT_, (gl_int_t)sh.memory);

    gl_int_t prev = 0;
    md_gl_fn.GetIntegerv(GL_TEXTURE_BINDING_2D_, &prev);
    md_gl_fn.GenTextures(1, &out->name);
    md_gl_fn.BindTexture(GL_TEXTURE_2D_, out->name);
    md_gl_fn.TexParameteri(GL_TEXTURE_2D_, GL_TEXTURE_TILING_EXT_, GL_OPTIMAL_TILING_EXT_);
    md_gl_fn.TexStorageMem2DEXT(GL_TEXTURE_2D_, 1, internal, (gl_sizei_t)sh.width, (gl_sizei_t)sh.height, out->memory, 0);
    md_gl_texture_params(GL_TEXTURE_2D_);
    md_gl_fn.BindTexture(GL_TEXTURE_2D_, (gl_uint_t)prev);
    out->target = GL_TEXTURE_2D_;

    if (sh.ready) {
        gl_uint_t sems[2] = {0, 0};
        md_gl_fn.GenSemaphoresEXT(2, sems);
        out->ready = sems[0];
        out->done  = sems[1];
        if (win32) {
            md_gl_fn.ImportSemaphoreWin32HandleEXT(out->ready, GL_HANDLE_TYPE_OPAQUE_WIN32_EXT_, (void*)(uintptr_t)sh.ready);
            md_gl_fn.ImportSemaphoreWin32HandleEXT(out->done,  GL_HANDLE_TYPE_OPAQUE_WIN32_EXT_, (void*)(uintptr_t)sh.done);
        } else {
            md_gl_fn.ImportSemaphoreFdEXT(out->ready, GL_HANDLE_TYPE_OPAQUE_FD_EXT_, (gl_int_t)sh.ready);
            md_gl_fn.ImportSemaphoreFdEXT(out->done,  GL_HANDLE_TYPE_OPAQUE_FD_EXT_, (gl_int_t)sh.done);
        }
    }

    const gl_enum_t err = md_gl_fn.GetError();
    if (err != GL_NO_ERROR_) {
        MD_LOG_ERROR("md_gpu_gl_import: GL error 0x%04X importing '%s'", err, desc->label ? desc->label : "texture");
        md_gpu_gl_release(out);
        return false;
    }
    return true;
}

void md_gpu_gl_begin_read(const md_gpu_gl_texture_t* t) {
    if (!t || !t->name) return;
    if (t->ready) {
        const gl_enum_t layout = GL_LAYOUT_GENERAL_EXT_;
        md_gl_fn.WaitSemaphoreEXT(t->ready, 0, NULL, 1, &t->name, &layout);
    }
    /* else: md_gpu's release returned once its writes were done */
}

void md_gpu_gl_end_read(const md_gpu_gl_texture_t* t) {
    if (!t || !t->name) return;
    if (t->done) {
        const gl_enum_t layout = GL_LAYOUT_GENERAL_EXT_;
        md_gl_fn.SignalSemaphoreEXT(t->done, 0, NULL, 1, &t->name, &layout);
        md_gl_fn.Flush();    /* the signal reaches the GPU before md_gpu's next acquire waits on it */
    } else {
        md_gl_fn.Finish();   /* GL cannot tell md_gpu: done on return */
    }
}

void md_gpu_gl_release(md_gpu_gl_texture_t* t) {
    if (!t) return;
    if (t->name) md_gl_fn.DeleteTextures(1, &t->name);
    if (t->memory && md_gl_fn.DeleteMemoryObjectsEXT) md_gl_fn.DeleteMemoryObjectsEXT(1, &t->memory);
    if (t->ready && md_gl_fn.DeleteSemaphoresEXT) {
        const gl_uint_t sems[2] = { t->ready, t->done };
        md_gl_fn.DeleteSemaphoresEXT(2, sems);
    }
    md_memset(t, 0, sizeof(*t));
}

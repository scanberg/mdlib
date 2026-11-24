#pragma once

/* md_gpu_gl: the OpenGL side of a texture md_gpu shares (MD_GPU_TEX_SHARED,
   see "Sharing a texture with OpenGL" in md_gpu.h): imported into the GL
   context current on the calling thread, and read there in turn with md_gpu's
   writes.

   No GL headers and no GL library: the entry points come from the caller's
   loader (glfwGetProcAddress, wglGetProcAddress, eglGetProcAddress, ...), so
   any program with a GL context can use it. Per platform it takes
   GL_EXT_memory_object with _fd or _win32 and GL_EXT_semaphore with _fd or
   _win32 (Vulkan), or CGLTexImageIOSurface2D (Metal; the loader finds the CGL
   entry points in OpenGL.framework).

   Each frame, after md_gpu_texture_share_release:

       md_gpu_gl_begin_read(&t);      // GL waits (on the GPU) for md_gpu's writes
       ... draw with t.name ...
       md_gpu_gl_end_read(&t);        // GL tells md_gpu it is done (and flushes)

   and only then md_gpu_texture_share_acquire for the next write. */

#include <core/md_gpu.h>

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef void* (*md_gpu_gl_loader_t)(const char* name);

typedef struct md_gpu_gl_texture_t {
    uint32_t name;              /* the GL texture                                     */
    uint32_t target;            /* GL_TEXTURE_2D, or GL_TEXTURE_RECTANGLE for an
                                   IOSurface (sampled in texels: sampler2DRect)       */
    uint32_t memory;            /* the GL memory object (Vulkan), else 0              */
    uint32_t ready, done;       /* the GL semaphores, else 0: then end_read waits on
                                   the host until GL is done (glFinish)               */
    md_gpu_share_kind_t kind;
    const void* fn;             /* the entry points, from the loader                  */
} md_gpu_gl_texture_t;

/* Imports `tex` into the current GL context. False, with the reason logged,
   when GL lacks what the kind needs, or runs on another GPU than the device
   (GL_DEVICE_UUID_EXT and GL_DRIVER_UUID_EXT against md_gpu's). Takes over the
   fds md_gpu_texture_share made for it: one import per share. */
bool md_gpu_gl_import(md_gpu_gl_texture_t* out, md_gpu_texture_t tex, md_gpu_gl_loader_t load);

void md_gpu_gl_begin_read(const md_gpu_gl_texture_t* t);
void md_gpu_gl_end_read(const md_gpu_gl_texture_t* t);

/* Deletes the GL texture, memory object and semaphores. Before or after the
   md_gpu texture is destroyed: they hold the memory of their own. */
void md_gpu_gl_release(md_gpu_gl_texture_t* t);

#ifdef __cplusplus
}
#endif

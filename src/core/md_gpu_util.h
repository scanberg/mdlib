/*
md_gpu_util.h

```
Utility layer built on top of md_gpu.h.

Provides higher-level helpers that are too opinionated for the core
abstraction but are broadly useful across different subsystems:

  - md_gpu_bump_t  : Persistent CPU-visible bump allocator for
                             staging uploads within a single command buffer.
```

Design notes:
  - Header-only; all functions are static inline.
  - No backend-specific code; depends only on md_gpu.h.
  - Thread-safety: none.  Each bump allocator instance is intended for
    use on a single thread (typically the render/compute thread).
*/

#pragma once

#include "md_gpu.h"

#include <core/md_common.h>

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* =============================
   Staging bump allocator
   =============================

   A grow-never-shrink CPU-visible buffer that is sub-divided by a
   monotonically advancing cursor.  Intended usage pattern per frame /
   per command recording session:

     1. md_gpu_bump_reset()          — cursor back to 0
     2. md_gpu_command_buffer_acquire()      — begin recording
     3. md_gpu_bump_alloc() × N     — reserve staging regions
        memcpy into alloc.cpu for each
        md_gpu_cmd_copy_buffer() for each   — record uploads
        md_gpu_cmd_barrier_buffer_ex()       — one barrier per target buf
     4. md_gpu_queue_submit()               — single submission

   Because each alloc occupies a distinct, non-overlapping region of the
   staging buffer the GPU can read all regions in the same submission
   without hazards.

   Overflow:
     md_gpu_bump_alloc() returns a zeroed struct on overflow.
     Call md_gpu_bump_grow() to expand the backing buffer; it
     preserves already-written data below the current cursor, so it is
     safe to call between alloc calls within the same recording session
     if you know you will not overflow again.  The simpler pattern is to
     pre-grow at session start once you know the total upload size.
*/

/* Result of a single sub-allocation.  cpu is a persistent CPU pointer
   valid until the staging buffer is destroyed or md_gpu_bump_grow()
   is called.  offset is the byte offset within the backing md_gpu_buffer_t
   to pass to md_gpu_cmd_copy_buffer() as src_offset. */
typedef struct md_gpu_alloc_t {
    void*  cpu;     /* CPU-writable pointer into the staging buffer */
    size_t offset;  /* byte offset within the backing buffer        */
} md_gpu_alloc_t;

typedef struct md_gpu_bump_alloc_t {
    md_gpu_device_t device;    /* backing device */
    md_gpu_buffer_t buffer;    /* backing CPU-visible buffer (owned) */
    uint8_t*        cpu_base;  /* == md_gpu_buffer_cpu_ptr(buffer)   */
    size_t          cursor;    /* next allocation starts here        */
    size_t          capacity;  /* total size of backing buffer       */
} md_gpu_bump_alloc_t;

/* Initialize or grow the backing buffer to at least `min_capacity` bytes.
   Safe to call on a zero-initialised struct (first-time init).
   Safe to call mid-session: existing data below cursor is preserved via
   memcpy into the new buffer.
   Returns false if the allocation fails; the struct is left unchanged. */
static inline bool md_gpu_bump_ensure(md_gpu_bump_alloc_t* bump, size_t min_capacity) {
    ASSERT(bump->device);
    if (bump->capacity >= min_capacity) return true;

    /* Round up to next power-of-two for headroom, minimum 64 KiB. */
    size_t new_cap = 65536;
    while (new_cap < min_capacity) new_cap *= 2;

    md_gpu_buffer_desc_t desc = { .size = new_cap, .flags = MD_GPU_BUFFER_CPU_VISIBLE };
    md_gpu_buffer_t new_buf = md_gpu_buffer_create(bump->device, &desc);
    if (!new_buf) return false;

    uint8_t* new_cpu = (uint8_t*)md_gpu_buffer_cpu_ptr(new_buf);

    /* Preserve already-written staging data so in-progress uploads remain valid. */
    if (bump->cursor > 0 && bump->cpu_base) {
        MEMCPY(new_cpu, bump->cpu_base, bump->cursor);
    }

    if (bump->buffer) md_gpu_buffer_destroy(bump->buffer);

    bump->buffer   = new_buf;
    bump->cpu_base = new_cpu;
    bump->capacity = new_cap;
    return true;
}

static inline bool md_gpu_bump_init(md_gpu_bump_alloc_t* bump, md_gpu_device_t device, size_t initial_capacity) {
    bump->device = device;
    return md_gpu_bump_ensure(bump, initial_capacity);
}

/* Free the backing buffer.  The struct is left zero-initialised and can
   be reused with md_gpu_bump_grow(). */
static inline void md_gpu_bump_free(md_gpu_bump_alloc_t* bump) {
    if (bump->buffer) md_gpu_buffer_destroy(bump->buffer);
    bump->device   = NULL;
    bump->buffer   = NULL;
    bump->cpu_base = NULL;
    bump->cursor   = 0;
    bump->capacity = 0;
}

/* Reset the cursor to 0.  Call once before starting a new command
   recording session.  Does not free or zero the backing memory. */
static inline void md_gpu_bump_reset(md_gpu_bump_alloc_t* bump) {
    bump->cursor = 0;
}

/* Convenience: grow if needed, then allocate.
   Combines md_gpu_bump_grow + md_gpu_bump_alloc.
   Returns a zeroed md_gpu_alloc_t on allocation failure. */
static inline md_gpu_alloc_t md_gpu_bump_push(md_gpu_bump_alloc_t* bump, size_t size) {
    const size_t align = 256; /* for GPU copy optimality */
    size_t off = (bump->cursor + align - 1) & ~(align - 1);
    if (off + size > bump->capacity) {
        if (!md_gpu_bump_ensure(bump, off + size)) {
            md_gpu_alloc_t empty = {0};
            return empty;
        }
    }
    bump->cursor = off + size;
    md_gpu_alloc_t result = { .cpu = bump->cpu_base + off, .offset = off };
    return result;
}

/* Total bytes consumed since the last reset. */
static inline size_t md_gpu_bump_used(const md_gpu_bump_alloc_t* bump) {
    return bump->cursor;
}

/* Remaining bytes available without a grow. */
static inline size_t md_gpu_bump_avail(const md_gpu_bump_alloc_t* bump) {
    return bump->capacity > bump->cursor ? bump->capacity - bump->cursor : 0;
}

#ifdef __cplusplus
}
#endif

#pragma once

// XDR (RFC 4506) primitives, as used by the GROMACS file formats (xtc, trr, edr, tpr).
//
// XDR is big endian and 4 byte aligned: every item occupies a multiple of 4 bytes, and variable
// length data (strings, opaque blobs) is padded with zeros up to the next multiple of 4.
//
// This is deliberately header only and allocation free. The formats that use it read a frame or a
// file into memory first and then decode it, so the only thing needed is a way to pull values out
// of a buffer, in two flavours:
//
//   md_xdr_load_*   Raw loads from a pointer. No bounds checking, for decoders that have already
//                   validated the size of what they are reading (the xtc frame decoder).
//   md_xdr_read_*   Reads through a cursor (md_xdr_t) that checks bounds. A read that would run past
//                   the end fails, leaves the cursor where it was, and marks the cursor as failed.
//                   The failure is sticky: every later read fails as well, so a parser can issue a
//                   sequence of reads and check md_xdr_ok() once at the end rather than after each.
//                   Values of failed reads are zeroed, never left uninitialized.

#include <core/md_common.h>
#include <core/md_str.h>

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// ### RAW LOADS ###

static inline uint32_t md_xdr_load_u32(const void* ptr) {
    uint32_t v;
    MEMCPY(&v, ptr, sizeof(v));
#if __LITTLE_ENDIAN__
    v = BSWAP32(v);
#endif
    return v;
}

static inline uint64_t md_xdr_load_u64(const void* ptr) {
    uint64_t v;
    MEMCPY(&v, ptr, sizeof(v));
#if __LITTLE_ENDIAN__
    v = BSWAP64(v);
#endif
    return v;
}

static inline int32_t md_xdr_load_i32(const void* ptr) {
    const uint32_t u = md_xdr_load_u32(ptr);
    int32_t v;
    MEMCPY(&v, &u, sizeof(v));
    return v;
}

static inline int64_t md_xdr_load_i64(const void* ptr) {
    const uint64_t u = md_xdr_load_u64(ptr);
    int64_t v;
    MEMCPY(&v, &u, sizeof(v));
    return v;
}

static inline float md_xdr_load_f32(const void* ptr) {
    const uint32_t u = md_xdr_load_u32(ptr);
    float v;
    MEMCPY(&v, &u, sizeof(v));
    return v;
}

static inline double md_xdr_load_f64(const void* ptr) {
    const uint64_t u = md_xdr_load_u64(ptr);
    double v;
    MEMCPY(&v, &u, sizeof(v));
    return v;
}

static inline void md_xdr_load_i32_array(int32_t* out, const void* ptr, size_t count) {
    const uint8_t* p = (const uint8_t*)ptr;
    for (size_t i = 0; i < count; ++i) {
        out[i] = md_xdr_load_i32(p + i * 4);
    }
}

static inline void md_xdr_load_f32_array(float* out, const void* ptr, size_t count) {
    const uint8_t* p = (const uint8_t*)ptr;
    for (size_t i = 0; i < count; ++i) {
        out[i] = md_xdr_load_f32(p + i * 4);
    }
}

static inline void md_xdr_load_f64_array(double* out, const void* ptr, size_t count) {
    const uint8_t* p = (const uint8_t*)ptr;
    for (size_t i = 0; i < count; ++i) {
        out[i] = md_xdr_load_f64(p + i * 8);
    }
}

// Bytes a variable length item of len bytes occupies, including its padding
static inline size_t md_xdr_padded_size(size_t len) {
    return (len + 3) & ~(size_t)3;
}

// ### CURSOR ###

typedef struct md_xdr_t {
    const uint8_t* data;
    size_t size;
    size_t pos;
    bool   error;
} md_xdr_t;

static inline md_xdr_t md_xdr_init(const void* data, size_t size) {
    md_xdr_t xdr = { (const uint8_t*)data, data ? size : 0, 0, false };
    return xdr;
}

static inline bool   md_xdr_ok(const md_xdr_t* xdr)        { return !xdr->error; }
static inline size_t md_xdr_remaining(const md_xdr_t* xdr) { return xdr->size - xdr->pos; }
static inline const uint8_t* md_xdr_ptr(const md_xdr_t* xdr) { return xdr->data + xdr->pos; }

// Reserves bytes at the cursor. Returns the pointer to them and advances, or NULL (and fails the
// cursor) if they are not there. This is what every read below is built on.
static inline const uint8_t* md_xdr_take(md_xdr_t* xdr, size_t bytes) {
    if (xdr->error || bytes > xdr->size - xdr->pos) {
        xdr->error = true;
        return NULL;
    }
    const uint8_t* ptr = xdr->data + xdr->pos;
    xdr->pos += bytes;
    return ptr;
}

static inline bool md_xdr_skip(md_xdr_t* xdr, size_t bytes) {
    return md_xdr_take(xdr, bytes) != NULL;
}

static inline bool md_xdr_read_u32(md_xdr_t* xdr, uint32_t* out) {
    const uint8_t* p = md_xdr_take(xdr, 4);
    *out = p ? md_xdr_load_u32(p) : 0;
    return p != NULL;
}

static inline bool md_xdr_read_i32(md_xdr_t* xdr, int32_t* out) {
    const uint8_t* p = md_xdr_take(xdr, 4);
    *out = p ? md_xdr_load_i32(p) : 0;
    return p != NULL;
}

// XDR hyper: 8 bytes, most significant word first
static inline bool md_xdr_read_u64(md_xdr_t* xdr, uint64_t* out) {
    const uint8_t* p = md_xdr_take(xdr, 8);
    *out = p ? md_xdr_load_u64(p) : 0;
    return p != NULL;
}

static inline bool md_xdr_read_i64(md_xdr_t* xdr, int64_t* out) {
    const uint8_t* p = md_xdr_take(xdr, 8);
    *out = p ? md_xdr_load_i64(p) : 0;
    return p != NULL;
}

static inline bool md_xdr_read_f32(md_xdr_t* xdr, float* out) {
    const uint8_t* p = md_xdr_take(xdr, 4);
    *out = p ? md_xdr_load_f32(p) : 0.0f;
    return p != NULL;
}

static inline bool md_xdr_read_f64(md_xdr_t* xdr, double* out) {
    const uint8_t* p = md_xdr_take(xdr, 8);
    *out = p ? md_xdr_load_f64(p) : 0.0;
    return p != NULL;
}

static inline bool md_xdr_read_i32_array(md_xdr_t* xdr, int32_t* out, size_t count) {
    const uint8_t* p = (count <= SIZE_MAX / 4) ? md_xdr_take(xdr, count * 4) : md_xdr_take(xdr, SIZE_MAX);
    if (!p) {
        MEMSET(out, 0, count * sizeof(int32_t));
        return false;
    }
    md_xdr_load_i32_array(out, p, count);
    return true;
}

static inline bool md_xdr_read_f32_array(md_xdr_t* xdr, float* out, size_t count) {
    const uint8_t* p = (count <= SIZE_MAX / 4) ? md_xdr_take(xdr, count * 4) : md_xdr_take(xdr, SIZE_MAX);
    if (!p) {
        MEMSET(out, 0, count * sizeof(float));
        return false;
    }
    md_xdr_load_f32_array(out, p, count);
    return true;
}

static inline bool md_xdr_read_f64_array(md_xdr_t* xdr, double* out, size_t count) {
    const uint8_t* p = (count <= SIZE_MAX / 8) ? md_xdr_take(xdr, count * 8) : md_xdr_take(xdr, SIZE_MAX);
    if (!p) {
        MEMSET(out, 0, count * sizeof(double));
        return false;
    }
    md_xdr_load_f64_array(out, p, count);
    return true;
}

// Fixed length opaque data: len bytes followed by padding to a multiple of 4. The bytes are not
// copied, *out points into the buffer.
static inline bool md_xdr_read_opaque(md_xdr_t* xdr, const uint8_t** out, size_t len) {
    const uint8_t* p = (len <= SIZE_MAX - 3) ? md_xdr_take(xdr, md_xdr_padded_size(len)) : md_xdr_take(xdr, SIZE_MAX);
    *out = p;
    return p != NULL;
}

// Variable length opaque data / string: a u32 byte count, the bytes, then padding. The result is a
// view into the buffer and is NOT zero terminated. max_len guards against a corrupt count asking
// for more than a sane string could hold; pass SIZE_MAX to accept anything that fits in the buffer.
static inline bool md_xdr_read_string(md_xdr_t* xdr, str_t* out, size_t max_len) {
    const size_t pos = xdr->pos;
    uint32_t len = 0;
    const uint8_t* p = NULL;
    if (md_xdr_read_u32(xdr, &len)) {
        if (len > max_len) {
            xdr->error = true;
        } else {
            md_xdr_read_opaque(xdr, &p, len);
        }
    }
    if (!p) {
        xdr->pos = pos;
        *out = (str_t){0};
        return false;
    }
    *out = (str_t){ (const char*)p, len };
    return true;
}

#ifdef __cplusplus
}
#endif

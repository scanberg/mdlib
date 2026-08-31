#pragma once

#include <stddef.h>

typedef enum {
    MD_COORD_STREAM_LAYOUT_SOA = 0,
    MD_COORD_STREAM_LAYOUT_AOS = 1,
} md_coord_stream_layout_t;

typedef struct md_coord_stream_t {
    md_coord_stream_layout_t layout;

    // Number of elements available in the stream (for asserts; optional but useful)
    size_t count;

    // Optional indirection (maps logical i -> source index)
    const int* idx;

    union {
        // SOA is always assumed to be tightly packed with no padding, so the stride is implicitly sizeof(float) and the offsets are implicitly 0.
        struct {
            const float* x;
            const float* y;
            const float* z;
        } soa;

        struct {
            // AOS is assumed to contain packed XYZ with some stride to the next triplet, i.e. a typical vertex struct.
            const float* base;
            size_t stride;   // in bytes
        } aos;
    };
} md_coord_stream_t;

static inline md_coord_stream_t md_coord_stream_from_soa(const float* x, const float* y, const float* z, const int* idx, size_t count) {
    md_coord_stream_t stream = {0};
    stream.layout = MD_COORD_STREAM_LAYOUT_SOA;
    stream.count = count;
    stream.idx = idx;
    stream.soa.x = x;
    stream.soa.y = y;
    stream.soa.z = z;
    return stream;
}

static inline md_coord_stream_t md_coord_stream_from_aos(const float* base_xyz, size_t stride_xyz, const int* idx, size_t count) {
    md_coord_stream_t stream = {0};
    stream.layout = MD_COORD_STREAM_LAYOUT_AOS;
    stream.count = count;
    stream.idx = idx;
    stream.aos.base = base_xyz;
    stream.aos.stride = stride_xyz;
    return stream;
}

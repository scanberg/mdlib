#include <md_hdf5_direct.h>

#include <md_attributes.h>
#include <md_hdf5.h>

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_simd.h>

#include <zlib.h>

#include <inttypes.h>
#include <string.h>

enum {
    DIRECT_NONE = 0,
    DIRECT_F32,
    DIRECT_F64,
    DIRECT_I8,
    DIRECT_U8,
    DIRECT_I16,
    DIRECT_U16,
    DIRECT_I32,
    DIRECT_U32,
    DIRECT_I64,
    DIRECT_U64,
};

// ## Describing the storage, with the lock held

static uint32_t direct_type(hid_t type, uint32_t* elem_size, uint32_t* big_endian) {
    *elem_size = 0;
    *big_endian = 0;
    if (type < 0) return DIRECT_NONE;
    const H5T_class_t cls = H5Tget_class(type);
    const size_t size = H5Tget_size(type);
    const H5T_order_t order = H5Tget_order(type);
    if (order != H5T_ORDER_LE && order != H5T_ORDER_BE) return DIRECT_NONE;
    uint32_t t = DIRECT_NONE;
    if (cls == H5T_FLOAT) {
        if      (H5Tequal(type, H5T_IEEE_F32LE) > 0 || H5Tequal(type, H5T_IEEE_F32BE) > 0) t = DIRECT_F32;
        else if (H5Tequal(type, H5T_IEEE_F64LE) > 0 || H5Tequal(type, H5T_IEEE_F64BE) > 0) t = DIRECT_F64;
    } else if (cls == H5T_INTEGER && H5Tget_precision(type) == size * 8 && H5Tget_offset(type) == 0) {
        const bool sign = H5Tget_sign(type) == H5T_SGN_2;
        switch (size) {
        case 1: t = sign ? DIRECT_I8  : DIRECT_U8;  break;
        case 2: t = sign ? DIRECT_I16 : DIRECT_U16; break;
        case 4: t = sign ? DIRECT_I32 : DIRECT_U32; break;
        case 8: t = sign ? DIRECT_I64 : DIRECT_U64; break;
        default: break;
        }
    }
    if (t != DIRECT_NONE) {
        *elem_size = (uint32_t)size;
        *big_endian = order == H5T_ORDER_BE;
    }
    return t;
}

static bool direct_filter_known(H5Z_filter_t f) {
    return f == H5Z_FILTER_DEFLATE || f == H5Z_FILTER_SHUFFLE || f == H5Z_FILTER_FLETCHER32;
}

// Where a chunk of the grid is in the table, from the offset of its first element
static uint64_t direct_chunk_index(const md_hdf5_direct_t* d, const hsize_t* offset) {
    uint64_t idx = 0;
    for (uint32_t i = 0; i < d->rank; ++i) idx = idx * d->grid[i] + offset[i] / d->chunk[i];
    return idx;
}

typedef struct direct_walk_t {
    md_hdf5_direct_t* d;
    hsize_t userblock;
    bool ok;
} direct_walk_t;

static void direct_walk_add(direct_walk_t* w, const hsize_t* offset, unsigned filter_mask, haddr_t addr, hsize_t size) {
    md_hdf5_direct_t* d = w->d;
    for (uint32_t i = 0; i < d->rank; ++i) {
        if (offset[i] % d->chunk[i] != 0 || offset[i] >= d->grid[i] * d->chunk[i]) { w->ok = false; return; }
    }
    if (addr == HADDR_UNDEF || size > UINT32_MAX) { w->ok = false; return; }
    const uint64_t idx = direct_chunk_index(d, offset);
    d->table[idx].offset = (int64_t)(w->userblock + addr);
    d->table[idx].size = (uint32_t)size;
    d->table[idx].filter_mask = filter_mask;
}

#if H5_VERSION_GE(1, 14, 0)
static int direct_chunk_visit(const hsize_t* offset, unsigned filter_mask, haddr_t addr, hsize_t size, void* data) {
    direct_walk_t* w = (direct_walk_t*)data;
    direct_walk_add(w, offset, filter_mask, addr, size);
    return w->ok ? H5_ITER_CONT : H5_ITER_STOP;
}
#endif

md_hdf5_direct_t* md_hdf5_direct_describe(hid_t file, hid_t dset, md_allocator_i* alloc) {
    MD_ASSERT(alloc);
    const hid_t space = H5Dget_space(dset);
    if (space < 0) return NULL;
    const int rank = H5Sget_simple_extent_ndims(space);
    hsize_t dims[MD_HDF5_DIRECT_MAX_RANK] = {0};
    const bool rank_ok = rank >= 1 && rank <= MD_HDF5_DIRECT_MAX_RANK && H5Sget_simple_extent_dims(space, dims, NULL) == rank;
    H5Sclose(space);
    if (!rank_ok) return NULL;

    md_hdf5_direct_t head = {0};
    head.rank = (uint32_t)rank;
    for (int i = 0; i < rank; ++i) head.dims[i] = dims[i];
    head.mode = MD_HDF5_DIRECT_HDF5;

    const hid_t type = H5Dget_type(dset);
    head.type = direct_type(type, &head.elem_size, &head.big_endian);
    if (type >= 0) H5Tclose(type);

    const hid_t dcpl = H5Dget_create_plist(dset);
    hsize_t userblock = 0;
    {
        const hid_t fcpl = H5Fget_create_plist(file);
        if (fcpl >= 0) {
            H5Pget_userblock(fcpl, &userblock);
            H5Pclose(fcpl);
        }
    }

    md_hdf5_direct_t* d = NULL;
    const H5D_layout_t layout = (dcpl >= 0 && head.type != DIRECT_NONE && H5Pget_external_count(dcpl) == 0) ? H5Pget_layout(dcpl) : H5D_LAYOUT_ERROR;
    if (layout == H5D_CONTIGUOUS) {
        const haddr_t addr = H5Dget_offset(dset);
        if (addr != HADDR_UNDEF) {
            head.mode = MD_HDF5_DIRECT_CONTIGUOUS;
            head.base = userblock + addr;
        }
    } else if (layout == H5D_CHUNKED) {
        hsize_t chunk[MD_HDF5_DIRECT_MAX_RANK] = {0};
        bool ok = H5Pget_chunk(dcpl, rank, chunk) == rank;
        const int num_filters = ok ? H5Pget_nfilters(dcpl) : -1;
        ok = ok && num_filters >= 0 && num_filters <= MD_HDF5_DIRECT_MAX_FILTERS;
        for (int i = 0; ok && i < num_filters; ++i) {
            unsigned int flags = 0, values[8];
            size_t num_values = MD_COUNTOF(values);
            char name[64];
            const H5Z_filter_t f = H5Pget_filter2(dcpl, (unsigned)i, &flags, &num_values, values, sizeof(name), name, NULL);
            ok = direct_filter_known(f);
            head.filter[i] = (uint32_t)f;
        }
        uint64_t num_chunks = 1;
        for (int i = 0; ok && i < rank; ++i) {
            ok = chunk[i] > 0;
            if (!ok) break;
            head.chunk[i] = chunk[i];
            head.grid[i] = (dims[i] + chunk[i] - 1) / chunk[i];
            num_chunks *= head.grid[i];
            ok = num_chunks < (1ull << 31);
        }
        if (ok) {
            head.num_filters = (uint32_t)num_filters;
            head.mode = MD_HDF5_DIRECT_CHUNKED;
            head.num_chunks = num_chunks;
            head.size = sizeof(md_hdf5_direct_t) + num_chunks * sizeof(md_hdf5_direct_chunk_t);
            d = md_alloc(alloc, head.size);
            if (!d) goto fail;
            *d = head;
            for (uint64_t c = 0; c < num_chunks; ++c) d->table[c] = (md_hdf5_direct_chunk_t){ .offset = -1 };

            direct_walk_t walk = { d, userblock, true };
#if H5_VERSION_GE(1, 14, 0)
            // One walk over the chunk index: asking for each chunk on its own costs a lookup each
            walk.ok = num_chunks == 0 || H5Dchunk_iter(dset, H5P_DEFAULT, direct_chunk_visit, &walk) >= 0;
#else
            // Before 1.14, chunk by chunk: 1.10 takes no H5S_ALL for the space
            const hid_t fspace = H5Dget_space(dset);
            hsize_t written = 0;
            walk.ok = fspace >= 0 && H5Dget_num_chunks(dset, fspace, &written) >= 0;
            for (hsize_t c = 0; walk.ok && c < written; ++c) {
                hsize_t offset[MD_HDF5_DIRECT_MAX_RANK] = {0};
                unsigned mask = 0;
                haddr_t addr = HADDR_UNDEF;
                hsize_t size = 0;
                walk.ok = H5Dget_chunk_info(dset, fspace, c, offset, &mask, &addr, &size) >= 0;
                if (walk.ok) direct_walk_add(&walk, offset, mask, addr, size);
            }
            if (fspace >= 0) H5Sclose(fspace);
#endif
            if (!walk.ok) {
                // HDF5 knows what to do with it
                md_free(alloc, d, head.size);
                d = NULL;
                head.mode = MD_HDF5_DIRECT_HDF5;
                head.num_chunks = 0;
                head.num_filters = 0;
            }
        }
    }
    if (!d) {
        if (head.mode != MD_HDF5_DIRECT_CONTIGUOUS) {
            head.mode = MD_HDF5_DIRECT_HDF5;
            head.num_chunks = 0;
            head.num_filters = 0;
        }
        head.size = sizeof(md_hdf5_direct_t);
        d = md_alloc(alloc, head.size);
        if (d) *d = head;
    }
fail:
    if (dcpl >= 0) H5Pclose(dcpl);
    return d;
}

// ## Reading

static size_t direct_row_length(const md_hdf5_direct_t* d) {
    size_t n = 1;
    for (uint32_t i = 1; i < d->rank; ++i) n *= (size_t)d->dims[i];
    return n;
}

static inline uint16_t direct_load16(const uint8_t* p, bool swap) { uint16_t u; MD_MEMCPY(&u, p, 2); return swap ? MD_BSWAP16(u) : u; }
static inline uint32_t direct_load32(const uint8_t* p, bool swap) { uint32_t u; MD_MEMCPY(&u, p, 4); return swap ? MD_BSWAP32(u) : u; }
static inline uint64_t direct_load64(const uint8_t* p, bool swap) { uint64_t u; MD_MEMCPY(&u, p, 8); return swap ? MD_BSWAP64(u) : u; }
static inline float    direct_f32(uint32_t u) { float f; MD_MEMCPY(&f, &u, 4); return f; }
static inline double   direct_f64(uint64_t u) { double f; MD_MEMCPY(&f, &u, 8); return f; }

// One loop per type: what a frame costs after it is inflated is mostly this
#define DIRECT_CONVERT(SIZE, EXPR) \
    for (size_t i = 0; i < count; ++i) { \
        const uint8_t* p = raw + i * (SIZE); \
        const double v = (double)(EXPR); \
        out[i] = (float)(v * scale); \
    }

static void direct_convert(float* out, const uint8_t* raw, size_t count, const md_hdf5_direct_t* d, double scale) {
    const bool swap = d->big_endian != 0;
    switch (d->type) {
    case DIRECT_F32:
        if (!swap && scale == 1.0) {
            MD_MEMCPY(out, raw, count * sizeof(float));
        } else if (!swap) {
            for (size_t i = 0; i < count; ++i) { float f; MD_MEMCPY(&f, raw + i * 4, 4); out[i] = (float)(f * scale); }
        } else {
            DIRECT_CONVERT(4, direct_f32(direct_load32(p, true)));
        }
        break;
    case DIRECT_F64: DIRECT_CONVERT(8, direct_f64(direct_load64(p, swap))); break;
    case DIRECT_I8:  DIRECT_CONVERT(1, (int8_t)p[0]); break;
    case DIRECT_U8:  DIRECT_CONVERT(1, p[0]); break;
    case DIRECT_I16: DIRECT_CONVERT(2, (int16_t)direct_load16(p, swap)); break;
    case DIRECT_U16: DIRECT_CONVERT(2, direct_load16(p, swap)); break;
    case DIRECT_I32: DIRECT_CONVERT(4, (int32_t)direct_load32(p, swap)); break;
    case DIRECT_U32: DIRECT_CONVERT(4, direct_load32(p, swap)); break;
    case DIRECT_I64: DIRECT_CONVERT(8, (int64_t)direct_load64(p, swap)); break;
    case DIRECT_U64: DIRECT_CONVERT(8, direct_load64(p, swap)); break;
    default: MD_MEMSET(out, 0, count * sizeof(float)); break;
    }
}

#undef DIRECT_CONVERT

static void direct_unshuffle(uint8_t* dst, const uint8_t* src, size_t bytes, size_t elem) {
    if (elem <= 1 || bytes < elem) {
        MD_MEMCPY(dst, src, bytes);
        return;
    }
    const size_t n = bytes / elem;
    // Element by element, so that the writes run through dst once; the reads are elem streams through src
    if (elem == 4) {
        const uint8_t *s0 = src, *s1 = src + n, *s2 = src + 2 * n, *s3 = src + 3 * n;
        for (size_t i = 0; i < n; ++i) {
            uint8_t* o = dst + i * 4;
            o[0] = s0[i]; o[1] = s1[i]; o[2] = s2[i]; o[3] = s3[i];
        }
    } else if (elem == 8) {
        for (size_t i = 0; i < n; ++i) {
            uint8_t* o = dst + i * 8;
            for (size_t b = 0; b < 8; ++b) o[b] = src[b * n + i];
        }
    } else {
        for (size_t b = 0; b < elem; ++b) {
            const uint8_t* s = src + b * n;
            for (size_t i = 0; i < n; ++i) dst[i * elem + b] = s[i];
        }
    }
    // What does not fill an element is stored as it is
    MD_MEMCPY(dst + n * elem, src + n * elem, bytes - n * elem);
}

static bool direct_inflate(uint8_t* dst, size_t* dst_len, size_t cap, const uint8_t* src, size_t src_len) {
    // A chunk is at most 4 GB, which is what zlib's counters hold
    if (src_len > UINT32_MAX || cap > UINT32_MAX) return false;
    z_stream z;
    MD_MEMSET(&z, 0, sizeof(z));
    if (inflateInit(&z) != Z_OK) return false;
    z.next_in   = (Bytef*)src;
    z.avail_in  = (uInt)src_len;
    z.next_out  = dst;
    z.avail_out = (uInt)cap;
    const int ret = inflate(&z, Z_FINISH);
    *dst_len = (size_t)z.total_out;
    inflateEnd(&z);
    return ret == Z_STREAM_END;
}

// ## Fletcher-32, as HDF5 computes it
//
// HDF5's (H5_checksum_fletcher32): the bytes as 16-bit big-endian words w_1.. w_n, a last odd byte the high half of a
// word, both sums from 0 and kept mod 65535 by folding the carry back in:
//
//     s1 = sum w_k,   s2 = sum (n - k + 1) w_k,   checksum = s2 << 16 | s1
//
// with one quirk of the folding: a sum that is a nonzero multiple of 65535 comes out as 0xffff, not 0.
//
// HDF5 computes it word by word, s2 += s1, a chain no compiler vectorizes. Here every lane of a vector keeps sums of
// its own - over rounds t of L lanes, word k = tL + j, lane j has v1_j = sum w and v2_j = the sum of v1_j before each
// round - and s1 = sum v1_j, s2 = L sum v2_j + sum (L - j) v1_j. Additions only, 16 words at a time.

#define FLETCHER_LANES  16
#define FLETCHER_ROUNDS 256     // per block: v2 stays below 2^32 (256 * 255 / 2 * 65535)

typedef struct direct_fletcher_t {
    uint64_t s1, s2;            // mod 65535
    bool     nonzero;           // any word nonzero: what tells 0 from 0xffff at the end
} direct_fletcher_t;

static inline void direct_fletcher_word(direct_fletcher_t* f, uint32_t w) {
    f->s1 += w;
    f->s2 += f->s1;
    f->nonzero |= w != 0;
}

uint32_t md_hdf5_fletcher32(const void* data, size_t bytes) {
    const uint8_t* p = (const uint8_t*)data;
    const size_t n = bytes / 2;     // whole words
    direct_fletcher_t f = {0};
    size_t k = 0;

    const size_t block = FLETCHER_LANES * FLETCHER_ROUNDS;
    md_256i any = md_mm256_setzero_si256();
    for (; k + block <= n; k += block) {
        md_256i v1a = md_mm256_setzero_si256(), v1b = v1a, v2a = v1a, v2b = v1a;
        const uint8_t* q = p + 2 * k;
        for (size_t t = 0; t < FLETCHER_ROUNDS; ++t, q += 2 * FLETCHER_LANES) {
            md_128i wa = md_mm_loadu_si128((const md_128i*)q);
            md_128i wb = md_mm_loadu_si128((const md_128i*)(q + 16));
            // Big endian words to the machine's order, widened to 32 bits: words 0-7 and 8-15
            wa = md_mm_or_si128(simde_mm_slli_epi16(wa, 8), simde_mm_srli_epi16(wa, 8));
            wb = md_mm_or_si128(simde_mm_slli_epi16(wb, 8), simde_mm_srli_epi16(wb, 8));
            const md_256i a = simde_mm256_cvtepu16_epi32(wa);
            const md_256i b = simde_mm256_cvtepu16_epi32(wb);
            v2a = md_mm256_add_epi32(v2a, v1a);
            v2b = md_mm256_add_epi32(v2b, v1b);
            v1a = md_mm256_add_epi32(v1a, a);
            v1b = md_mm256_add_epi32(v1b, b);
            any = md_mm256_or_si256(any, md_mm256_or_si256(a, b));
        }
        uint32_t v1[FLETCHER_LANES], v2[FLETCHER_LANES];
        md_mm256_storeu_si256((md_256i*)v1, v1a);
        md_mm256_storeu_si256((md_256i*)(v1 + 8), v1b);
        md_mm256_storeu_si256((md_256i*)v2, v2a);
        md_mm256_storeu_si256((md_256i*)(v2 + 8), v2b);
        uint64_t a = 0, c = 0;
        for (size_t j = 0; j < FLETCHER_LANES; ++j) {
            a += v1[j];
            c += (uint64_t)FLETCHER_LANES * v2[j] + (uint64_t)(FLETCHER_LANES - j) * v1[j];
        }
        // The block after what came before it: each of its words also adds s1 so far to s2
        f.s2 = (f.s2 + (uint64_t)block * f.s1 + c) % 65535;
        f.s1 = (f.s1 + a) % 65535;
    }
    {
        uint32_t lanes[8];
        md_mm256_storeu_si256((md_256i*)lanes, any);
        for (size_t j = 0; j < 8; ++j) f.nonzero |= lanes[j] != 0;
    }
    for (; k < n; ++k) {
        direct_fletcher_word(&f, ((uint32_t)p[2 * k] << 8) | p[2 * k + 1]);
        if ((k & 0xfff) == 0xfff) { f.s1 %= 65535; f.s2 %= 65535; }
    }
    if (bytes & 1) direct_fletcher_word(&f, (uint32_t)p[bytes - 1] << 8);

    uint32_t s1 = (uint32_t)(f.s1 % 65535), s2 = (uint32_t)(f.s2 % 65535);
    if (f.nonzero && s1 == 0) s1 = 0xffff;
    if (f.nonzero && s2 == 0) s2 = 0xffff;
    return (s2 << 16) | s1;
}

// Whether the four bytes after data hold its checksum: as HDF5 writes it (a little endian uint32), or as HDF5 before
// 1.6.3 did on little endian machines, each half's two bytes swapped - which HDF5 still accepts
static bool direct_fletcher_check(const uint8_t* data, size_t bytes) {
    const uint8_t* c = data + bytes;
    const uint32_t stored = (uint32_t)c[0] | ((uint32_t)c[1] << 8) | ((uint32_t)c[2] << 16) | ((uint32_t)c[3] << 24);
    const uint32_t sum = md_hdf5_fletcher32(data, bytes);
    const uint32_t swapped = ((sum & 0x00ff00ffu) << 8) | ((sum >> 8) & 0x00ff00ffu);
    return stored == sum || stored == swapped;
}

// The bytes of a chunk as HDF5 would hand them out: the filters it was written through undone, in reverse. The last
// step writes into out, which holds out_bytes and not a byte more. *bad_checksum is set when a checksum does not match.
static bool direct_decode(uint8_t* out, size_t out_bytes, const uint8_t* stored, size_t stored_bytes, const md_hdf5_direct_t* d, uint32_t mask, bool* bad_checksum, md_temp_scope_t temp) {
    // The last filter to undo that writes anew (inflate, unshuffle) writes into out, the ones before it into scratch;
    // unless a checksum is still to be taken off after it, and out would be four bytes short
    int final = -1;
    for (int i = 0; i < (int)d->num_filters; ++i) {
        if (mask & (1u << i)) continue;
        if (d->filter[i] == H5Z_FILTER_FLETCHER32) break;
        final = i;
        break;
    }
    const size_t cap = out_bytes + 64;
    uint8_t* buf[2] = { NULL, NULL };
    const uint8_t* cur = stored;
    size_t len = stored_bytes;
    int next = 0;
    for (int i = (int)d->num_filters - 1; i >= 0; --i) {
        if (mask & (1u << i)) continue;
        uint8_t* dst = out;
        size_t dst_cap = out_bytes;
        if (i != final && d->filter[i] != H5Z_FILTER_FLETCHER32) {
            if (!buf[next]) buf[next] = md_temp_alloc(temp, cap);
            if (!buf[next]) return false;
            dst = buf[next];
            dst_cap = cap;
            next ^= 1;
        }
        switch (d->filter[i]) {
        case H5Z_FILTER_FLETCHER32:
            // The checksum is the last four bytes, of the bytes before it, which were just read or just decoded
            if (len < 4) return false;
            len -= 4;
            if (!direct_fletcher_check(cur, len)) {
                *bad_checksum = true;
                return false;
            }
            break;
        case H5Z_FILTER_DEFLATE: {
            size_t n = 0;
            if (!direct_inflate(dst, &n, dst_cap, cur, len)) return false;
            cur = dst; len = n;
            break;
        }
        case H5Z_FILTER_SHUFFLE:
            if (len > dst_cap) return false;
            direct_unshuffle(dst, cur, len, d->elem_size);
            cur = dst;
            break;
        default:
            return false;
        }
    }
    if (len != out_bytes) return false;
    if (cur != out) MD_MEMCPY(out, cur, len);
    return true;
}

// The decoded bytes of chunk c, kept by io when there is one
static const uint8_t* direct_chunk(const md_hdf5_direct_t* d, uint64_t c, size_t chunk_bytes, md_str_t path, md_attribute_io_t* io, md_temp_scope_t temp) {
    const md_hdf5_direct_chunk_t* e = &d->table[c];
    const uint64_t key = md_attribute_io_block_key(path, e->offset);
    bool hit = false;
    uint8_t* out = md_attribute_io_block(io, key, chunk_bytes, &hit);
    if (out && hit) return out;
    const bool kept = out != NULL;
    if (!out) out = md_temp_alloc(temp, chunk_bytes);
    if (!out) return NULL;

    bool active = false;
    for (uint32_t i = 0; i < d->num_filters; ++i) active |= !(e->filter_mask & (1u << i));
    bool ok;
    if (!active) {
        ok = e->size == chunk_bytes && md_attribute_io_read_at(io, path, e->offset, out, chunk_bytes) == chunk_bytes;
    } else {
        uint8_t* stored = md_temp_alloc(temp, e->size);
        bool bad_checksum = false;
        ok = stored && md_attribute_io_read_at(io, path, e->offset, stored, e->size) == e->size &&
             direct_decode(out, chunk_bytes, stored, e->size, d, e->filter_mask, &bad_checksum, temp);
        if (bad_checksum) {
            MD_LOG_ERROR("HDF5: the chunk at %" PRId64 " of '" MD_STR_FMT "' does not match its Fletcher-32 checksum: the file is damaged", e->offset, MD_STR_ARG(path));
        }
    }
    if (!ok) {
        if (kept) md_attribute_io_block_forget(io, key);
        return NULL;
    }
    return out;
}

// Through HDF5, under the lock: the rows of axis 1 the range touches, then the range of them
static bool direct_read_hdf5(float* out, const md_hdf5_direct_t* d, md_str_t path, const char* dataset, size_t frame, size_t first, size_t count, double scale) {
    size_t inner = 1;
    for (uint32_t i = 2; i < d->rank; ++i) inner *= (size_t)d->dims[i];
    hsize_t start[MD_HDF5_DIRECT_MAX_RANK] = {0}, extent[MD_HDF5_DIRECT_MAX_RANK] = {0};
    start[0] = frame;
    extent[0] = 1;
    size_t skip = 0, total = 1;
    if (d->rank >= 2) {
        const size_t row0 = first / inner;
        const size_t row1 = (first + count - 1) / inner;
        start[1] = row0;
        extent[1] = row1 - row0 + 1;
        skip = first - row0 * inner;
        total = (size_t)extent[1] * inner;
        for (uint32_t i = 2; i < d->rank; ++i) extent[i] = d->dims[i];
    }

    char cpath[4096];
    md_str_copy_to_char_buf(cpath, sizeof(cpath), path);
    md_temp_scope_t temp = md_temp_begin();
    float* buf = total == count ? out : md_temp_alloc_array(temp, float, total);
    bool ok = false;
    if (buf) {
        const md_hdf5_lock_t lock = md_hdf5_lock();
        const hid_t file = H5Fopen(cpath, H5F_ACC_RDONLY, H5P_DEFAULT);
        const hid_t dset = file >= 0 ? H5Dopen2(file, dataset, H5P_DEFAULT) : -1;
        if (dset >= 0) {
            const hid_t space = H5Dget_space(dset);
            const hsize_t n = total;
            const hid_t mem = H5Screate_simple(1, &n, NULL);
            ok = H5Sselect_hyperslab(space, H5S_SELECT_SET, start, NULL, extent, NULL) >= 0 &&
                 H5Dread(dset, H5T_NATIVE_FLOAT, mem, space, H5P_DEFAULT, buf) >= 0;
            H5Sclose(mem);
            H5Sclose(space);
            H5Dclose(dset);
        }
        if (file >= 0) H5Fclose(file);
        md_hdf5_unlock(lock);
    }
    if (ok) {
        if (buf != out) MD_MEMCPY(out, buf + skip, count * sizeof(float));
        if (scale != 1.0) for (size_t i = 0; i < count; ++i) out[i] = (float)(out[i] * scale);
    }
    md_temp_end(temp);
    return ok;
}

// Steps k through [lo[i], hi[i]) over the axes [beg, end), the last fastest; false when it has been through them all
static bool direct_next(size_t* k, const size_t* lo, const size_t* hi, uint32_t beg, uint32_t end) {
    for (uint32_t i = end; i-- > beg;) {
        if (++k[i] < hi[i]) return true;
        k[i] = lo[i];
    }
    return false;
}

// The chunked read: every chunk of the frame's slab whose elements reach into the range, and of each the runs along
// the last axis that do. 1 when read, 0 when it failed, -1 when a chunk it needs was never written.
static int direct_read_chunked(float* out, const md_hdf5_direct_t* d, md_str_t path, size_t frame, size_t first, size_t count, double scale, md_attribute_io_t* io, md_temp_scope_t temp) {
    const uint32_t r = d->rank;
    const uint32_t L = r - 1;           // the last axis
    const size_t last = first + count;  // exclusive
    // Strides of the frame's row, flattened over axes 1.. L
    size_t stride[MD_HDF5_DIRECT_MAX_RANK] = {0};
    {
        size_t s = 1;
        for (uint32_t i = L; i >= 1; --i) { stride[i] = s; s *= (size_t)d->dims[i]; }
    }
    size_t chunk_elems = 1;
    for (uint32_t i = 0; i < r; ++i) chunk_elems *= (size_t)d->chunk[i];
    const size_t chunk_bytes = chunk_elems * d->elem_size;
    const size_t g0 = frame / (size_t)d->chunk[0];
    const size_t k0 = frame - g0 * (size_t)d->chunk[0];

    size_t zero[MD_HDF5_DIRECT_MAX_RANK] = {0}, grid[MD_HDF5_DIRECT_MAX_RANK] = {0};
    for (uint32_t i = 1; i < r; ++i) grid[i] = (size_t)d->grid[i];

    size_t g[MD_HDF5_DIRECT_MAX_RANK] = {0};    // the chunk, over axes 1.. L
    do {
        // Its extent in the dataset, and the first and last element of it in the row
        size_t cs[MD_HDF5_DIRECT_MAX_RANK] = {0}, ce[MD_HDF5_DIRECT_MAX_RANK] = {0};
        size_t lo = 0, hi = 0;
        for (uint32_t i = 1; i < r; ++i) {
            cs[i] = g[i] * (size_t)d->chunk[i];
            ce[i] = MD_MIN(cs[i] + (size_t)d->chunk[i], (size_t)d->dims[i]);
            lo += cs[i] * stride[i];
            hi += (ce[i] - 1) * stride[i];
        }
        if (hi < first || lo >= last) continue;

        uint64_t c = g0;
        for (uint32_t i = 1; i < r; ++i) c = c * d->grid[i] + g[i];
        if (d->table[c].offset < 0) return -1;
        const uint8_t* bytes = direct_chunk(d, c, chunk_bytes, path, io, temp);
        if (!bytes) return 0;

        if (r == 1) {
            direct_convert(out, bytes + k0 * d->elem_size, 1, d, scale);
            continue;
        }
        // The runs along the last axis, one for each element of the chunk's axes 1.. L-1
        const size_t run = ce[L] - cs[L];
        size_t k[MD_HDF5_DIRECT_MAX_RANK] = {0};
        for (uint32_t i = 1; i < L; ++i) k[i] = cs[i];
        do {
            size_t t = cs[L] * stride[L];   // where the run is in the row
            size_t local = k0;              // and in the chunk
            for (uint32_t i = 1; i < L; ++i) {
                t += k[i] * stride[i];
                local = local * (size_t)d->chunk[i] + (k[i] - cs[i]);
            }
            local *= (size_t)d->chunk[L];
            const size_t a = MD_MAX(t, first), b = MD_MIN(t + run, last);
            if (a < b) direct_convert(out + (a - first), bytes + (local + (a - t)) * d->elem_size, b - a, d, scale);
        } while (direct_next(k, cs, ce, 1, L));
    } while (direct_next(g, zero, grid, 1, r));
    return 1;
}

bool md_hdf5_direct_read_f32(float* out, const md_hdf5_direct_t* d, md_str_t path, const char* dataset, size_t frame, size_t first, size_t count, double scale, md_attribute_io_t* io) {
    if (!d || !out || frame >= d->dims[0] || count == 0) return false;
    const size_t row = direct_row_length(d);
    if (first + count > row) return false;

    if (d->mode == MD_HDF5_DIRECT_CONTIGUOUS) {
        md_temp_scope_t temp = md_temp_begin();
        const size_t bytes = count * d->elem_size;
        uint8_t* raw = md_temp_alloc(temp, bytes);
        const int64_t at = (int64_t)(d->base + ((uint64_t)frame * row + first) * d->elem_size);
        const bool ok = raw && md_attribute_io_read_at(io, path, at, raw, bytes) == bytes;
        if (ok) direct_convert(out, raw, count, d, scale);
        md_temp_end(temp);
        return ok;
    }
    if (d->mode == MD_HDF5_DIRECT_CHUNKED) {
        md_temp_scope_t temp = md_temp_begin();
        const int r = direct_read_chunked(out, d, path, frame, first, count, scale, io, temp);
        md_temp_end(temp);
        if (r >= 0) return r == 1;
    }
    return dataset && direct_read_hdf5(out, d, path, dataset, frame, first, count, scale);
}

void md_hdf5_direct_frame_spans(int64_t* offsets, int64_t* sizes, const md_hdf5_direct_t* d, size_t num_frames) {
    const size_t row_bytes = d ? direct_row_length(d) * d->elem_size : 0;
    bool whole = d && d->mode == MD_HDF5_DIRECT_CHUNKED;
    for (uint32_t i = 1; whole && i < d->rank; ++i) whole = d->grid[i] == 1;
    for (size_t f = 0; f < num_frames; ++f) {
        offsets[f] = -1;
        sizes[f] = 0;
        if (!d || f >= d->dims[0]) continue;
        if (d->mode == MD_HDF5_DIRECT_CONTIGUOUS) {
            offsets[f] = (int64_t)(d->base + f * row_bytes);
            sizes[f] = (int64_t)row_bytes;
        } else if (whole) {
            const md_hdf5_direct_chunk_t* e = &d->table[f / d->chunk[0]];
            bool filtered = false;
            for (uint32_t i = 0; i < d->num_filters; ++i) filtered |= !(e->filter_mask & (1u << i));
            if (e->offset < 0) {
                continue;
            } else if (!filtered && e->size >= d->chunk[0] * row_bytes) {
                // Frames stored as they are, one after the other
                offsets[f] = e->offset + (int64_t)((f % d->chunk[0]) * row_bytes);
                sizes[f] = (int64_t)row_bytes;
            } else {
                // The chunk, which holds the frame among others
                offsets[f] = e->offset;
                sizes[f] = e->size;
            }
        }
    }
}

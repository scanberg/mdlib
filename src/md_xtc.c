#include <md_xtc.h>
#include <md_xdr.h>

#include <md_system.h>

#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_vec_math.h>

#include <libdivide.h>
#include <stdio.h>

#if defined(__SIZEOF_INT128__)
#define HAS_INT128_T
#endif

#define MD_XTC_CACHE_MAGIC   0x8281237612371
#define MD_XTC_CACHE_VERSION 6   // 6: the shared run cache header; 5: the step and box of every frame

#define XTC_MAGIC 1995

/* XTC small header size (natoms<=9).
*  > int(4) magic
*  > int(4) natoms
*  > int(4) step
*  > float(4) time
*  > 9xfloat(4) box
*  > int(4) natoms (again)
*/
#define XTC_SMALL_HEADER_SIZE 56

/* Size of uncompressed coordinates for one atom.
* 3xfloat(4) x
*/
#define XTC_SMALL_COORDS_SIZE 12

/* XTC header size (natoms>=10).
* Compressed trajectories contain some additional values:
*  > float(4) precision
*  > 3xint(4) minint
*  > 3xint(4) maxint
*  > int(4) smallidx
* See `xdrfile_compress_coord_double()`.
*/
#define XTC_HEADER_SIZE (XTC_SMALL_HEADER_SIZE + 32)

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(magicints) / sizeof(*magicints))

#define DIV_T struct libdivide_u64_branchfree_t
#define DIV_INIT(x) libdivide_u64_branchfree_gen(x)
#define DIV(x, y) libdivide_u64_branchfree_do(x, y)

typedef md_128i v4i_t;
#define v4i_set(x, y, z, w) md_mm_set_epi32(w, z, y, x)
#define v4i_set1(x)         md_mm_set1_epi32(x)
#define v4i_add(a, b)       md_mm_add_epi32(a, b)
#define v4i_sub(a, b)       md_mm_sub_epi32(a, b)
#define v4i_load(addr)      md_mm_loadu_epi32(addr)

// Number of guard bytes that must be readable past the end of a bitstream so
// the stateless 16-byte unaligned reads in extract_bits_be_raw_* never overrun.
#define MD_XTC_STREAM_GUARD_BYTES 16
// Required alignment (in bytes) for the bitstream buffer.
#define MD_XTC_STREAM_ALIGNMENT   16

typedef struct bit_data_t {
    uint32_t num_of_bits;
    uint32_t part_mask;
    uint32_t big_shift;
    uint32_t sml_shift;
} bit_data_t;

typedef struct unpack_data_t {
    uint32_t   size_y;
    uint32_t   size_z;
    bit_data_t bit;
    DIV_T      div_zy;
    DIV_T      div_z;
} unpack_data_t;

static const uint32_t magicints[] = {
    0,        0,        0,       0,       0,       0,       0,       0,       0,       8,
    10,       12,       16,      20,      25,      32,      40,      50,      64,      80,
    101,      128,      161,     203,     256,     322,     406,     512,     645,     812,
    1024,     1290,     1625,    2048,    2580,    3250,    4096,    5060,    6501,    8192,
    10321,    13003,    16384,   20642,   26007,   32768,   41285,   52015,   65536,   82570,
    104031,   131072,   165140,  208063,  262144,  330280,  416127,  524287,  660561,  832255,
    1048576,  1321122,  1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042, 8388607,
    10568983, 13316085, 16777216};

// These contain precalculated multiplies and shifts for the corresponding magic numbers from FIRSTIDX to LASTIDX
static const DIV_T denoms_64_1[] = {
    {0, 2}, {0X999999999999999A, 3}, {0X5555555555555556, 3}, {0, 3}, {0X999999999999999A, 4}, {0X47AE147AE147AE15, 4}, {0, 4}, {0X999999999999999A, 5}, {0X47AE147AE147AE15, 5}, {0, 5},
    {0X999999999999999A, 6}, {0X446F86562D9FAEE5, 6}, {0, 6}, {0X970E4F80CB8727C1, 7}, {0X42D6625D51F86EFA, 7}, {0, 7}, {0X970E4F80CB8727C1, 8}, {0X42D6625D51F86EFA, 8}, {0, 8},
    {0X966CC01966CC0197, 9}, {0X42D6625D51F86EFA, 9}, {0, 9}, {0X966CC01966CC0197, 10}, {0X42A386615BD84CE1, 10}, {0, 10}, {0X966CC01966CC0197, 11}, {0X42A386615BD84CE1, 11}, {0, 11},
    {0X9E74F8832302A17F, 12}, {0X4296D1E340EF6BF0, 12}, {0, 12}, {0X9662AB64ED3938ED, 13}, {0X4290780440A7950F, 13}, {0, 13}, {0X9662AB64ED3938ED, 14}, {0X428D4B2CC2273C78, 14}, {0, 14},
    {0X9660264BCF9BE4F1, 15}, {0X428BB4C7031B065B, 15}, {0, 15}, {0X9660264BCF9BE4F1, 16}, {0X428AE995A39B4B45, 16}, {0, 16}, {0X9660264BCF9BE4F1, 17}, {0X428A83FD53DC320A, 17}, {0, 17},
    {0X9660264BCF9BE4F1, 18}, {0X428A513143FCBC0B, 18}, {0X200004000081, 18}, {0X965FFDFA81C292FF, 19}, {0X428A37CB420D0364, 19}, {0, 19}, {0X965FFDFA81C292FF, 20},
    {0X428A37CB420D0364, 20}, {0, 20}, {0X965FF3E62F8C3EC9, 21}, {0X428A3171C2311550, 21}, {0, 21}, {0X965FEEDC06A114B2, 22}, {0X428A3171C2311550, 22}, {0X20000040001, 22},
    {0X965FEC56F2377FA7, 23}, {0X428A2FDB624419CC, 23}, {0, 23}, 
};

// Precalculated multiply and shifts for squared magic numbers
static const DIV_T denoms_64_2[] = {
    {0, 5}, {0X47AE147AE147AE15, 6}, {0XC71C71C71C71C71D, 7}, {0, 7}, {0X47AE147AE147AE15, 8}, {0XA36E2EB1C432CA58, 9}, {0, 9}, {0X47AE147AE147AE15, 10}, {0XA36E2EB1C432CA58, 11},
    {0, 11}, {0X47AE147AE147AE15, 12}, {0X9B2A7C9FE8B617F1, 13}, {0, 13}, {0X439F40CC28F760CC, 14}, {0X972002FB5C05974D, 15}, {0, 15}, {0X439F40CC28F760CC, 16}, {0X972002FB5C05974D, 17},
    {0, 17}, {0X429E8FC19BD1F6F1, 18}, {0X972002FB5C05974D, 19}, {0, 19}, {0X429E8FC19BD1F6F1, 20}, {0X969FC68151912805, 21}, {0, 21}, {0X429E8FC19BD1F6F1, 22}, {0X969FC68151912805, 23},
    {0, 23}, {0X4F7F449D26A949DD, 24}, {0X967FC0DA51BF3DDD, 25}, {0, 25}, {0X428E8ED5E63B8E8A, 26}, {0X966FBF71EFB24C2F, 27}, {0, 27}, {0X428E8ED5E63B8E8A, 28}, {0X9667BF187DCBD4ED, 29},
    {0, 29}, {0X428A8ECA9A8C6397, 30}, {0X9663BF027395D269, 31}, {0, 31}, {0X428A8ECA9A8C6397, 32}, {0X9661BEFD1A08C80F, 33}, {0, 33}, {0X428A8ECA9A8C6397, 34}, {0X9660BEFBD82195C9, 35},
    {0, 35}, {0X428A8ECA9A8C6397, 36}, {0X96603EFB91E54C07, 37}, {0X40000C000201, 37}, {0X428A4ECA87BD4923, 38}, {0X965FFEFB8574EA53, 39}, {0, 39}, {0X428A4ECA87BD4923, 40},
    {0X965FFEFB8574EA53, 41}, {0, 41}, {0X428A3ECA8603773D, 42}, {0X965FEEFB84B5956C, 43}, {0, 43}, {0X428A36CA8598D950, 44}, {0X965FEEFB84B5956C, 45}, {0X400000C0001, 45},
    {0X428A32CA85801D1A, 46}, {0X965FEAFB84AB8C66, 47}, {0, 47},
};

// =====================================================================
// Stateless bit extraction helpers.
//
// The decoder consumes a big-endian XTC bitstream sequentially. Reads are
// performed directly from the underlying byte buffer using the current bit
// offset; no sliding window / cached qwords are kept. This avoids a long
// dependency chain through a stateful bitreader and exposes much more
// instruction-level parallelism (each extract is independent given an offset).
//
// All extractors load up to 16 bytes from `base + (bit_offset >> 3)`. The
// caller MUST ensure that `MD_XTC_STREAM_GUARD_BYTES` (16) bytes past the
// logical end of the bitstream are readable.
// =====================================================================

static FORCE_INLINE int sizeofint(int size) {
    unsigned int num = 1;
    int num_of_bits = 0;

    while (size >= (int)num && num_of_bits < 32) {
        num_of_bits++;
        num *= 2;
    }
    return num_of_bits;
}

static FORCE_INLINE int sizeofints(int num_of_ints, unsigned int sizes[]) {
    unsigned int num_of_bytes, num_of_bits, bytes[32], bytecnt, tmp, num;
    num_of_bytes = 1;
    bytes[0] = 1;
    num_of_bits = 0;
    for (int i = 0; i < num_of_ints; i++) {
        tmp = 0;
        for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
            tmp = bytes[bytecnt] * sizes[i] + tmp;
            bytes[bytecnt] = tmp & 0xff;
            tmp >>= 8;
        }
        while (tmp != 0) {
            bytes[bytecnt++] = tmp & 0xff;
            tmp >>= 8;
        }
        num_of_bytes = bytecnt;
    }
    num = 1;
    num_of_bytes--;
    while (bytes[num_of_bytes] >= num) {
        num_of_bits++;
        num *= 2;
    }
    return num_of_bits + num_of_bytes * 8;
}

static FORCE_INLINE void write_coord(float* dst, v4i_t coord, md_128 invp) {
    md_128 data = md_mm_mul_ps(md_mm_cvtepi32_ps(coord), invp);
    MEMCPY(dst, &data, 3 * sizeof(float));
}

static FORCE_INLINE void write_coord_soa(float* RESTRICT x, float* RESTRICT y, float* RESTRICT z, size_t idx, v4i_t coord, md_128 scale) {
    ALIGNAS(16) float data[4];
    md_mm_store_ps(data, md_mm_mul_ps(md_mm_cvtepi32_ps(coord), scale));
    x[idx] = data[0];
    y[idx] = data[1];
    z[idx] = data[2];
}

static FORCE_INLINE void init_unpack_bit_data(bit_data_t* data, uint32_t num_of_bits) {
    uint32_t partbits = num_of_bits & 7;
    data->num_of_bits = num_of_bits;
    data->big_shift   = 64 - ((num_of_bits + 7) & ~7);  // Align to next multiple of 8
    data->sml_shift   = (8 - partbits) & 7;
    data->part_mask   = partbits ? (1u << partbits) - 1u : 0xFFu;
}

// Fast path: extracts up to 57 bits via a single unaligned 8-byte read.
// Works because (bit_in_byte <= 7) + bit_length <= 64.
static FORCE_INLINE uint64_t extract_bits_be_raw_57(const uint8_t* base, size_t bit_offset, size_t bit_length) {
    size_t byte_offset = bit_offset >> 3;
    size_t bit_in_byte = bit_offset & 7;

    uint64_t raw;
    MEMCPY(&raw, base + byte_offset, sizeof(raw));
#if __LITTLE_ENDIAN__
    raw = BSWAP64(raw);
#endif
    return (raw << bit_in_byte) >> (64 - bit_length);
}

// Extracts 1..64 bits using a 16-byte unaligned read so that any bit
// alignment is handled correctly.
static FORCE_INLINE uint64_t extract_bits_be_raw_64(const uint8_t* base, size_t bit_offset, size_t bit_length) {
    size_t byte_offset = bit_offset >> 3;
    size_t bit_in_byte = bit_offset & 7;

    uint64_t raw[2];
    MEMCPY(raw, base + byte_offset, sizeof(raw));
#if __LITTLE_ENDIAN__
    raw[0] = BSWAP64(raw[0]);
    raw[1] = BSWAP64(raw[1]);
#endif
    uint64_t nz_mask  = -(uint64_t)(bit_in_byte != 0);
    uint64_t combined = (raw[0] << bit_in_byte) | ((raw[1] >> (64 - bit_in_byte)) & nz_mask);
    return combined >> (64 - bit_length);
}

// Extracts 65..121 bits. Output: out[0] = upper portion (still byte-shifted),
// out[1] = lower 64 bits (right-justified). Designed to feed unpack_coord128.
static FORCE_INLINE void extract_bits_be_raw_121(uint64_t out[2], const uint8_t* base, size_t bit_offset, size_t bit_length) {
    ASSERT(64 < bit_length && bit_length <= 121);

    size_t byte_offset = bit_offset >> 3;
    size_t bit_in_byte = bit_offset & 7;
    int    shift_right = (int)(128 - bit_length);

    uint64_t raw[2];
    MEMCPY(raw, base + byte_offset, sizeof(raw));
#if __LITTLE_ENDIAN__
    raw[0] = BSWAP64(raw[0]);
    raw[1] = BSWAP64(raw[1]);
#endif
    uint64_t nz_mask = -(uint64_t)(bit_in_byte != 0);
    uint64_t hi = (raw[0] << bit_in_byte);
    uint64_t lo = (raw[1] << bit_in_byte) | ((raw[0] >> (64 - bit_in_byte)) & nz_mask);
    out[0] = hi;
    out[1] = lo >> shift_right;
}

// Extracts 1..32 bits via a single unaligned 8-byte read.
static FORCE_INLINE uint32_t extract_bits_be_raw_32(const uint8_t* base, size_t bit_offset, size_t bit_length) {
    size_t byte_offset = bit_offset >> 3;
    size_t bit_in_byte = bit_offset & 7;

    uint64_t raw;
    MEMCPY(&raw, base + byte_offset, sizeof(raw));
#if __LITTLE_ENDIAN__
    raw = BSWAP64(raw);
#endif
    return (uint32_t)((raw << bit_in_byte) >> (64 - bit_length));
}

// Peek-style helper for the run-length 6-bit field. Reads up to 32 bits.
static FORCE_INLINE uint32_t extract_bits_be_raw_25(const uint8_t* base, size_t bit_offset, size_t bit_length) {
    size_t byte_offset = bit_offset >> 3;
    size_t bit_in_byte = bit_offset & 7;

    uint32_t raw;
    MEMCPY(&raw, base + byte_offset, sizeof(raw));
#if __LITTLE_ENDIAN__
    raw = BSWAP32(raw);
#endif
    return (raw << bit_in_byte) >> (32 - bit_length);
}

static FORCE_INLINE uint32_t unpack_uint32(uint32_t w, const bit_data_t* bits) {
    uint32_t t = ((w << bits->sml_shift) & (~0xFFu)) | (w & bits->part_mask);
    return BSWAP32(t) >> bits->big_shift;
}

// Combine a 64-bit packed word into (x, y, z) via two parallel divisions.
// The two libdivide ops are issued from the same source `v`, which shortens
// the critical dependency chain compared to a serial (x = v/zy, q = v - x*zy,
// y = q/z) formulation.
static FORCE_INLINE v4i_t unpack_coord64(uint64_t w, const unpack_data_t* unpack) {
    uint64_t t = ((w << unpack->bit.sml_shift) & (~0xFFull)) | (w & unpack->bit.part_mask);
    uint64_t v = BSWAP64(t) >> unpack->bit.big_shift;

    uint32_t x  = (uint32_t)DIV(v, &unpack->div_zy);
    uint64_t yz = (uint64_t)DIV(v, &unpack->div_z);

    uint32_t y = (uint32_t)(yz - (uint64_t)x  * unpack->size_y);
    uint32_t z = (uint32_t)(v  - (uint64_t)yz * unpack->size_z);

    return v4i_set(x, y, z, 0);
}

static FORCE_INLINE v4i_t unpack_coord128(uint64_t w[2], const unpack_data_t* unpack) {
    const uint64_t zy = (uint64_t)unpack->size_z * unpack->size_y;

    // w[0] is the still-shifted upper portion (see extract_bits_be_raw_121),
    // w[1] is the right-justified lower 64 bits but stored big-endian.
    uint64_t hi = ((w[1] << unpack->bit.sml_shift) & (~0xFFull)) | (w[1] & unpack->bit.part_mask);
    hi = BSWAP64(hi) >> unpack->bit.big_shift;
    uint64_t lo = BSWAP64(w[0]);

#ifdef HAS_INT128_T
    __uint128_t v = (__uint128_t)lo | ((__uint128_t)hi << 64);
    uint32_t x = (uint32_t)(v / zy);
    uint64_t q = (uint64_t)(v - x*zy);
    uint32_t y = (uint32_t)(q / unpack->size_z);
    uint32_t z = (uint32_t)(q % unpack->size_z);
    return v4i_set(x, y, z, 0);
#elif defined(_MSC_VER) && defined(_M_X64) && (_MSC_VER >= 1920)
    uint64_t q = 0;
    uint32_t x = (uint32_t)_udiv128(hi, lo, zy, &q);
    uint32_t y = (uint32_t)(q / unpack->size_z);
    uint32_t z = (uint32_t)(q % unpack->size_z);
    return v4i_set(x, y, z, 0);
#else
    // Generic schoolbook division across 32-bit limbs.
    const int      fullbytes = unpack->bit.num_of_bits >> 3;
    const int      partbits  = unpack->bit.num_of_bits  & 7;
    const int      num_of_bytes = fullbytes + (partbits ? 1 : 0);
    const uint32_t sizes[3] = {0, unpack->size_y, unpack->size_z};
    uint32_t nums[4]  = {0};
    uint32_t limbs[4] = {0};

    MEMCPY(limbs + 0, &lo, sizeof(uint32_t) * 2);
    MEMCPY(limbs + 2, &hi, sizeof(uint32_t) * 2);

    const int num_limbs = (num_of_bytes + 3) / 4;
    for (int i = 2; i > 0; --i) {
        const uint32_t d = sizes[i];
        uint32_t rem = 0;
        for (int j = num_limbs - 1; j >= 0; --j) {
            uint64_t cur = ((uint64_t)rem << 32) | (uint64_t)limbs[j];
            limbs[j] = (uint32_t)(cur / d);
            rem      = (uint32_t)(cur % d);
        }
        nums[i] = rem;
    }
    nums[0] = limbs[0];
    return v4i_load(nums);
#endif
}

// Picks the cheapest extraction kernel for `bit_length` and feeds unpack_coord*.
// `bit_length` may exceed 64 (rare). The <= 57 path is the hot one and is
// inlined first so the compiler can favor it.
static FORCE_INLINE v4i_t extract_and_unpack(const uint8_t* base, size_t bit_offset, size_t bit_length, const unpack_data_t* unpack) {
    if (bit_length <= 57) {
        uint64_t w = extract_bits_be_raw_57(base, bit_offset, bit_length);
        return unpack_coord64(w, unpack);
    } else if (bit_length <= 64) {
        uint64_t w = extract_bits_be_raw_64(base, bit_offset, bit_length);
        return unpack_coord64(w, unpack);
    } else {
        uint64_t w[2];
        extract_bits_be_raw_121(w, base, bit_offset, bit_length);
        return unpack_coord128(w, unpack);
    }
}

// Reads three independently bit-packed integers (one per axis). Used only
// when sizeints exceed the 24-bit combined limit.
static FORCE_INLINE v4i_t extract_ints3(const uint8_t* base, size_t bit_offset, const bit_data_t bits[3]) {
    uint32_t v[4] = {0};
    for (int i = 0; i < 3; ++i) {
        uint32_t num_bits = bits[i].num_of_bits;
        uint32_t w = extract_bits_be_raw_32(base, bit_offset, num_bits);
        v[i] = unpack_uint32(w, &bits[i]);
        bit_offset += num_bits;
    }
    return v4i_load(v);
}

static inline bool decode_header(const uint8_t* frame_ptr, md_xtc_header_t* out_header) {
    // Extract header
    int magic;
    magic = md_xdr_load_i32(frame_ptr);
    if (magic != XTC_MAGIC) {
        MD_LOG_ERROR("XTC: Magic number did not match");
        return false;
	}

    out_header->natoms = md_xdr_load_i32(frame_ptr + 4);
    out_header->step = md_xdr_load_i32(frame_ptr + 8);
	out_header->time = md_xdr_load_f32(frame_ptr + 12);
	md_xdr_load_f32_array((float*)out_header->box, frame_ptr + 16, 9);

    return true;
}

// The scan reads every frame header on its way through the file, and a header holds the step and
// the box as well as the time. steps and boxes may be NULL; boxes takes nine floats per frame, nm.
static size_t xtc_scan(md_file_t xdr, md_array(int64_t)* frame_offsets, md_array(double)* frame_times,
                       md_array(int64_t)* frame_steps, md_array(float)* frame_boxes, md_allocator_i* alloc);

size_t md_xtc_read_frame_offsets_and_times(md_file_t xdr, md_array(int64_t)* frame_offsets, md_array(double)* frame_times, md_allocator_i* alloc) {
    return xtc_scan(xdr, frame_offsets, frame_times, NULL, NULL, alloc);
}

static void xtc_scan_push_header(md_array(int64_t)* frame_steps, md_array(float)* frame_boxes, const md_xtc_header_t* h, md_allocator_i* alloc) {
    if (frame_steps) {
        md_array_push(*frame_steps, (int64_t)h->step, alloc);
    }
    if (frame_boxes) {
        const float* box = &h->box[0][0];
        for (int i = 0; i < 9; ++i) {
            md_array_push(*frame_boxes, box[i], alloc);
        }
    }
}

static size_t xtc_scan(md_file_t xdr, md_array(int64_t)* frame_offsets, md_array(double)* frame_times,
                       md_array(int64_t)* frame_steps, md_array(float)* frame_boxes, md_allocator_i* alloc) {
    size_t filesize = (size_t)md_file_size(xdr);

    if (filesize == 0) {
        MD_LOG_ERROR("XTC: Failed extract filesize");
        return 0;
    }

    /* Go to file beg */
    if (!md_file_seek(xdr, 0, MD_FILE_BEG)) {
        MD_LOG_ERROR("XTC: Failed to seek to beginning of file");
        return 0;
    }

    uint8_t frame_header_data[XTC_HEADER_SIZE];
    md_xtc_header_t xtc_header;
    size_t read_bytes;

    read_bytes = md_file_read(xdr, frame_header_data, XTC_SMALL_HEADER_SIZE);
    if (read_bytes != XTC_SMALL_HEADER_SIZE || !decode_header(frame_header_data, &xtc_header)) {
        MD_LOG_ERROR("XTC: File does not appear to be a valid xtc trajectory");
        return 0;
    }

    if (xtc_header.natoms <= 0) {
        MD_LOG_ERROR("XTC: Invalid number of atoms in header");
        return 0;
    }

    size_t num_frames = 0;

    // Push first frame
    md_array_push(*frame_offsets, 0, alloc);
    md_array_push(*frame_times, xtc_header.time, alloc);
    xtc_scan_push_header(frame_steps, frame_boxes, &xtc_header, alloc);
    num_frames += 1;

    /* Dont bother with compression for nine atoms or less */
    if (xtc_header.natoms <= 9) {
        const size_t framebytes = XTC_SMALL_HEADER_SIZE + XTC_SMALL_COORDS_SIZE * xtc_header.natoms;
        const size_t est_frames = (filesize / framebytes); /* Should we complain if framesize doesn't divide filesize? */

        md_array_ensure(*frame_offsets, est_frames, alloc);
        md_array_ensure(*frame_times,   est_frames, alloc);

        for (size_t i = 1; i < est_frames; i++) {
            const size_t offset = i * framebytes;
            
            bool success = false;
            if (md_file_seek(xdr, offset, MD_FILE_BEG)) {
                read_bytes = md_file_read(xdr, frame_header_data, XTC_SMALL_HEADER_SIZE);
                success = (read_bytes == XTC_SMALL_HEADER_SIZE) && decode_header(frame_header_data, &xtc_header);
            }

            // Push frame i
            if (success) {
                md_array_push(*frame_offsets, offset, alloc);
                md_array_push(*frame_times, xtc_header.time, alloc);
                xtc_scan_push_header(frame_steps, frame_boxes, &xtc_header, alloc);
                num_frames += 1;
            } else {
               MD_LOG_DEBUG("XTC: encountered corrupted frame header");
               break;
            }
        }
        md_array_push(*frame_offsets, num_frames * framebytes, alloc);
    } else {
        int framebytes = 0;
        int est_nframes = 0;

        /* Move pos back to end of first header */
        if (!md_file_seek(xdr, XTC_HEADER_SIZE, MD_FILE_BEG)) {
            return 0;
        }

        uint8_t framebytes_raw[4];
        if (md_file_read(xdr, framebytes_raw, sizeof(framebytes_raw)) != sizeof(framebytes_raw)) {
            MD_LOG_ERROR("XTC: Failed to read framebytes");
            return 0;
        }
        framebytes = (int)md_xdr_padded_size((size_t)md_xdr_load_u32(framebytes_raw)); /* Rounding to the next 32-bit boundary */

        /* Skip `framebytes` */
        if (!md_file_seek(xdr, framebytes, MD_FILE_CUR)) {
            MD_LOG_DEBUG("XTC: encountered corrupted frame");
            return 0;
        }

        est_nframes = (int)(filesize / (framebytes + XTC_HEADER_SIZE) + 1);
        /* First `framebytes` might be larger than average, so we would underestimate `est_nframes`*/
        est_nframes += est_nframes / 5;

        md_array_ensure(*frame_offsets, (size_t)est_nframes, alloc);
        md_array_ensure(*frame_times,   (size_t)est_nframes, alloc);

        while (true) {
            const int64_t offset = md_file_tell(xdr);
            if (offset == (int64_t)filesize) {
                // Good exit
                break;
            }

            read_bytes = md_file_read(xdr, frame_header_data, XTC_HEADER_SIZE);
            if (read_bytes != XTC_HEADER_SIZE || !decode_header(frame_header_data, &xtc_header)) {
                MD_LOG_DEBUG("XTC: encountered corrupted frame header");
                goto done;
            }

            /* Read how much to skip */
            if (md_file_read(xdr, framebytes_raw, sizeof(framebytes_raw)) != sizeof(framebytes_raw)) {
                MD_LOG_ERROR("XTC: Failed to read framebytes");
                goto done;
            }
            framebytes = (int)md_xdr_padded_size((size_t)md_xdr_load_u32(framebytes_raw)); /* Rounding to the next 32-bit boundary */

            /* Skip `framebytes` to next header */
            if (!md_file_seek(xdr, framebytes, MD_FILE_CUR)) {
                MD_LOG_DEBUG("XTC: encountered corrupted frame");
                goto done;
            }

            /* Store position in `offsets`, adjust for header */
            md_array_push(*frame_offsets, offset, alloc);
            md_array_push(*frame_times, xtc_header.time, alloc);
            xtc_scan_push_header(frame_steps, frame_boxes, &xtc_header, alloc);
            num_frames += 1;
        }
        // Add last offset
        md_array_push(*frame_offsets, filesize, alloc);
    }
done:
    return num_frames;
}

// xyz packed, the layout the file holds them in, scaled on the way out exactly as the SoA variant
// below scales: scale 10 turns the file's nm into Angstrom with no second pass over the coordinates.
static bool xtc_decode_frame_data_scaled(const uint8_t* frame_ptr, size_t frame_bytes, md_xtc_header_t* out_header, float* out_coords, size_t num_atoms, float scale) {
    if (frame_ptr == NULL || frame_bytes == 0) {
        return false;
    }

    if (frame_bytes < XTC_SMALL_HEADER_SIZE) {
        MD_LOG_ERROR("XTC: Frame size is too small to contain header");
        return false;
	}

    if (out_header) {
        if (!decode_header(frame_ptr, out_header)) {
            return false;
        }
        float* box = &out_header->box[0][0];
        for (int i = 0; i < 9; ++i) {
            box[i] *= scale;
        }
	}

    if (!out_coords) {
        return true;
    }

    int natoms;
	natoms = md_xdr_load_i32(frame_ptr + XTC_SMALL_HEADER_SIZE - 4);

    if (natoms != (int)num_atoms) {
        MD_LOG_ERROR("XTC: Number of atoms in frame header does not match expected number of atoms");
        return false;
    }

    size_t offset = XTC_SMALL_HEADER_SIZE;
    if (natoms <= 9) {
		// No compression for 9 atoms or less, just read the coordinates directly
		md_xdr_load_f32_array(out_coords, frame_ptr + offset, (size_t)natoms * 3);
        for (size_t i = 0; i < (size_t)natoms * 3; ++i) {
            out_coords[i] *= scale;
        }
        return true;
    }

    float precision;
    int32_t minint[3], maxint[3], smallidx;
	precision = md_xdr_load_f32(frame_ptr + offset); offset += 4;
	md_xdr_load_i32_array(minint, frame_ptr + offset, 3); offset += 12;
	md_xdr_load_i32_array(maxint, frame_ptr + offset, 3); offset += 12;
	smallidx = md_xdr_load_i32(frame_ptr + offset); offset += 4;

    uint32_t sizeint[3] = {
        (uint32_t)(maxint[0] - minint[0] + 1),
        (uint32_t)(maxint[1] - minint[1] + 1),
        (uint32_t)(maxint[2] - minint[2] + 1),
    };

    uint32_t bitsize = 0;
    bit_data_t bitsizeint[3] = {0};
    unpack_data_t big_unpack = {0};

    if ((sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffff) {
        init_unpack_bit_data(&bitsizeint[0], sizeofint(sizeint[0]));
        init_unpack_bit_data(&bitsizeint[1], sizeofint(sizeint[1]));
        init_unpack_bit_data(&bitsizeint[2], sizeofint(sizeint[2]));
    } else {
        bitsize = sizeofints(3, sizeint);
        big_unpack.size_y = sizeint[1];
        big_unpack.size_z = sizeint[2];
        big_unpack.div_zy = DIV_INIT((uint64_t)sizeint[1] * sizeint[2]);
        big_unpack.div_z  = DIV_INIT(sizeint[2]);
        init_unpack_bit_data(&big_unpack.bit, bitsize);
    }

    int idx = MAX(smallidx - 1, FIRSTIDX);
    int smaller = magicints[idx] / 2;
    int smallnum = magicints[smallidx] / 2;
    uint32_t smallsize = magicints[smallidx];

    unpack_data_t sml_unpack = {
        .size_y = smallsize,
        .size_z = smallsize,
        .div_zy = denoms_64_2[smallidx - FIRSTIDX],
        .div_z  = denoms_64_1[smallidx - FIRSTIDX],
    };
    init_unpack_bit_data(&sml_unpack.bit, smallidx);

    /* length in bytes */
    int32_t num_bytes = 0;
    num_bytes = md_xdr_load_i32(frame_ptr + offset); offset += 4;
    (void)num_bytes;

    const uint8_t* stream = frame_ptr + offset;
    size_t bit_offset = 0;

    float* lfp = out_coords;
    md_128 invp = md_mm_set1_ps(scale / precision);
    v4i_t vminint = v4i_set(minint[0], minint[1], minint[2], 0);
    v4i_t thiscoord;
    int run = 0;
    int run_count = 0;
    int atom_idx = 0;

    while (atom_idx < natoms) {
        if (bitsize == 0) {
            thiscoord = extract_ints3(stream, bit_offset, bitsizeint);
            bit_offset += (size_t)bitsizeint[0].num_of_bits + bitsizeint[1].num_of_bits + bitsizeint[2].num_of_bits;
        } else {
            thiscoord = extract_and_unpack(stream, bit_offset, bitsize, &big_unpack);
            bit_offset += bitsize;
        }

        thiscoord = v4i_add(thiscoord, vminint);

        uint32_t data = extract_bits_be_raw_25(stream, bit_offset, 6);
        uint32_t flag = data & 32;
        uint32_t skip = flag ? 6 : 1;
        bit_offset += skip;

        int is_smaller = 0;
        if (flag) {
            run = data & 31;
            run_count  = run / 3;
            is_smaller = run % 3;
            run -= is_smaller;
            is_smaller--;
        }

        int batch_size = run_count + 1;
        if (atom_idx + batch_size > natoms) {
            MD_LOG_ERROR("XTC: Buffer overrun during decompression.");
            goto done;
        }
        atom_idx += batch_size;

        if (run > 0) {
            v4i_t prevcoord = thiscoord;
            v4i_t vsmall = v4i_set1(smallnum);
            uint32_t sml_bits = sml_unpack.bit.num_of_bits;

            uint64_t w = extract_bits_be_raw_57(stream, bit_offset, sml_bits);
            v4i_t coord = unpack_coord64(w, &sml_unpack);
            bit_offset += sml_bits;
            thiscoord = v4i_add(coord, v4i_sub(thiscoord, vsmall));

            write_coord(lfp, thiscoord, invp); lfp += 3;
            write_coord(lfp, prevcoord, invp); lfp += 3;

            for (int i = 1; i < run_count; ++i) {
                w = extract_bits_be_raw_57(stream, bit_offset, sml_bits);
                coord = unpack_coord64(w, &sml_unpack);
                bit_offset += sml_bits;
                thiscoord = v4i_add(coord, v4i_sub(thiscoord, vsmall));
                write_coord(lfp, thiscoord, invp); lfp += 3;
            }
        } else {
            write_coord(lfp, thiscoord, invp); lfp += 3;
        }
        smallidx += is_smaller;
        if (is_smaller < 0) {
            smallnum = smaller;
            smaller = (smallidx > FIRSTIDX) ? magicints[smallidx - 1] / 2 : 0;
        } else if (is_smaller > 0) {
            smaller = smallnum;
            smallnum = magicints[smallidx] / 2;
        }
        if (smallidx < FIRSTIDX) {
            MD_LOG_ERROR("XTC: Invalid size found in 'xdrfile_decompress_coord_float'.");
            goto done;
        }
        if ((uint32_t)smallidx != sml_unpack.bit.num_of_bits) {
            uint32_t sml_size       = magicints[smallidx];
            sml_unpack.size_y       = sml_size;
            sml_unpack.size_z       = sml_size;
            sml_unpack.div_zy       = denoms_64_2[smallidx - FIRSTIDX];
            sml_unpack.div_z        = denoms_64_1[smallidx - FIRSTIDX];
            init_unpack_bit_data(&sml_unpack.bit, smallidx);
        }
    }

done:
    return atom_idx == natoms;
}

bool md_xtc_decode_frame_data(const uint8_t* frame_ptr, size_t frame_bytes, md_xtc_header_t* out_header, float* out_coords, size_t num_atoms) {
    return xtc_decode_frame_data_scaled(frame_ptr, frame_bytes, out_header, out_coords, num_atoms, 1.0f);
}

static bool md_xtc_decode_frame_data_soa_scaled(const uint8_t* frame_ptr, size_t frame_bytes, md_xtc_header_t* out_header, float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z, size_t num_atoms, float scale) {
    if (frame_ptr == NULL || frame_bytes == 0) {
        return false;
    }

    if (frame_bytes < XTC_SMALL_HEADER_SIZE) {
        MD_LOG_ERROR("XTC: Frame size is too small to contain header");
        return false;
    }

    md_xtc_header_t header = {0};
    if (!decode_header(frame_ptr, &header)) {
        return false;
    }

    if (out_header) {
        // Scale box dimensions with scale factor
        float* box = (float*)header.box;
        for (int i = 0; i < 9; ++i) {
            box[i] *= scale;
        }
        MEMCPY(out_header, &header, sizeof(md_xtc_header_t));
    }

    if (!out_x || !out_y || !out_z) {
        return true;
    }

    int natoms;
    natoms = md_xdr_load_i32(frame_ptr + XTC_SMALL_HEADER_SIZE - 4);

    if (natoms != (int)num_atoms) {
        MD_LOG_ERROR("XTC: Number of atoms in frame header does not match expected number of atoms");
        return false;
    }

    size_t offset = XTC_SMALL_HEADER_SIZE;
    if (natoms <= 9) {
        for (int i = 0; i < natoms; ++i) {
            float coord[3];
            md_xdr_load_f32_array(coord, frame_ptr + offset + (size_t)i * XTC_SMALL_COORDS_SIZE, 3);
            out_x[i] = coord[0] * scale;
            out_y[i] = coord[1] * scale;
            out_z[i] = coord[2] * scale;
        }
        return true;
    }

    float precision;
    int32_t minint[3], maxint[3], smallidx;
    precision = md_xdr_load_f32(frame_ptr + offset); offset += 4;
    md_xdr_load_i32_array(minint, frame_ptr + offset, 3); offset += 12;
    md_xdr_load_i32_array(maxint, frame_ptr + offset, 3); offset += 12;
    smallidx = md_xdr_load_i32(frame_ptr + offset); offset += 4;

    uint32_t sizeint[3] = {
        (uint32_t)(maxint[0] - minint[0] + 1),
        (uint32_t)(maxint[1] - minint[1] + 1),
        (uint32_t)(maxint[2] - minint[2] + 1),
    };

    uint32_t bitsize = 0;
    bit_data_t bitsizeint[3] = {0};
    unpack_data_t big_unpack = {0};

    if ((sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffff) {
        init_unpack_bit_data(&bitsizeint[0], sizeofint(sizeint[0]));
        init_unpack_bit_data(&bitsizeint[1], sizeofint(sizeint[1]));
        init_unpack_bit_data(&bitsizeint[2], sizeofint(sizeint[2]));
    } else {
        bitsize = sizeofints(3, sizeint);
        big_unpack.size_y = sizeint[1];
        big_unpack.size_z = sizeint[2];
        big_unpack.div_zy = DIV_INIT((uint64_t)sizeint[1] * sizeint[2]);
        big_unpack.div_z  = DIV_INIT(sizeint[2]);
        init_unpack_bit_data(&big_unpack.bit, bitsize);
    }

    int idx = MAX(smallidx - 1, FIRSTIDX);
    int smaller = magicints[idx] / 2;
    int smallnum = magicints[smallidx] / 2;
    uint32_t smallsize = magicints[smallidx];

    unpack_data_t sml_unpack = {
        .size_y = smallsize,
        .size_z = smallsize,
        .div_zy = denoms_64_2[smallidx - FIRSTIDX],
        .div_z  = denoms_64_1[smallidx - FIRSTIDX],
    };
    init_unpack_bit_data(&sml_unpack.bit, smallidx);

    int32_t num_bytes = 0;
    num_bytes = md_xdr_load_i32(frame_ptr + offset); offset += 4;
    (void)num_bytes;

    const uint8_t* stream = frame_ptr + offset;
    size_t bit_offset = 0;

    md_128 coord_scale = md_mm_set1_ps(scale / precision);
    v4i_t vminint = v4i_set(minint[0], minint[1], minint[2], 0);
    v4i_t thiscoord;
    int run = 0;
    int run_count = 0;
    int atom_idx = 0;

    while (atom_idx < natoms) {
        if (bitsize == 0) {
            thiscoord = extract_ints3(stream, bit_offset, bitsizeint);
            bit_offset += (size_t)bitsizeint[0].num_of_bits + bitsizeint[1].num_of_bits + bitsizeint[2].num_of_bits;
        } else {
            thiscoord = extract_and_unpack(stream, bit_offset, bitsize, &big_unpack);
            bit_offset += bitsize;
        }

        thiscoord = v4i_add(thiscoord, vminint);

        uint32_t data = extract_bits_be_raw_25(stream, bit_offset, 6);
        uint32_t flag = data & 32;
        uint32_t skip = flag ? 6 : 1;
        bit_offset += skip;

        int is_smaller = 0;
        if (flag) {
            run = data & 31;
            run_count  = run / 3;
            is_smaller = run % 3;
            run -= is_smaller;
            is_smaller--;
        }

        int batch_size = run_count + 1;
        if (atom_idx + batch_size > natoms) {
            MD_LOG_ERROR("XTC: Buffer overrun during decompression.");
            goto done;
        }

        if (run > 0) {
            v4i_t prevcoord = thiscoord;
            v4i_t vsmall = v4i_set1(smallnum);
            uint32_t sml_bits = sml_unpack.bit.num_of_bits;

            uint64_t w = extract_bits_be_raw_57(stream, bit_offset, sml_bits);
            v4i_t coord = unpack_coord64(w, &sml_unpack);
            bit_offset += sml_bits;
            thiscoord = v4i_add(coord, v4i_sub(thiscoord, vsmall));

            write_coord_soa(out_x, out_y, out_z, atom_idx++, thiscoord, coord_scale);
            write_coord_soa(out_x, out_y, out_z, atom_idx++, prevcoord, coord_scale);

            for (int i = 1; i < run_count; ++i) {
                w = extract_bits_be_raw_57(stream, bit_offset, sml_bits);
                coord = unpack_coord64(w, &sml_unpack);
                bit_offset += sml_bits;
                thiscoord = v4i_add(coord, v4i_sub(thiscoord, vsmall));
                write_coord_soa(out_x, out_y, out_z, atom_idx++, thiscoord, coord_scale);
            }
        } else {
            write_coord_soa(out_x, out_y, out_z, atom_idx++, thiscoord, coord_scale);
        }

        smallidx += is_smaller;
        if (is_smaller < 0) {
            smallnum = smaller;
            smaller = (smallidx > FIRSTIDX) ? magicints[smallidx - 1] / 2 : 0;
        } else if (is_smaller > 0) {
            smaller = smallnum;
            smallnum = magicints[smallidx] / 2;
        }
        if (smallidx < FIRSTIDX) {
            MD_LOG_ERROR("XTC: Invalid size found in 'xdrfile_decompress_coord_float'.");
            goto done;
        }
        if ((uint32_t)smallidx != sml_unpack.bit.num_of_bits) {
            uint32_t sml_size       = magicints[smallidx];
            sml_unpack.size_y       = sml_size;
            sml_unpack.size_z       = sml_size;
            sml_unpack.div_zy       = denoms_64_2[smallidx - FIRSTIDX];
            sml_unpack.div_z        = denoms_64_1[smallidx - FIRSTIDX];
            init_unpack_bit_data(&sml_unpack.bit, smallidx);
        }
    }

done:
    return atom_idx == natoms;
}

bool md_xtc_decode_frame_data_soa(const uint8_t* frame_ptr, size_t frame_bytes, md_xtc_header_t* out_header, float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z, size_t num_atoms) {
    return md_xtc_decode_frame_data_soa_scaled(frame_ptr, frame_bytes, out_header, out_x, out_y, out_z, num_atoms, 1.0f);
}

// Everything the scan learns about the file, which is everything short of the coordinates: where
// each frame is, and the time, step and box its header states. Kept beside the trajectory as
// '<file>.cache' so the scan is paid once per file rather than once per load.
typedef struct xtc_cache_t {
    md_run_cache_header_t header;
    int64_t* frame_offsets;   // num_frames + 1: the last one is the file size
    double*  frame_times;     // ps
    int64_t* frame_steps;
    float*   frame_boxes;     // 9 per frame, nm, row major
} xtc_cache_t;

static bool cache_read_block(md_file_t file, void** dst, size_t bytes, md_allocator_i* alloc) {
    *dst = md_alloc(alloc, bytes);
    return *dst && md_file_read(file, *dst, bytes) == bytes;
}

// On failure nothing is freed here: alloc is always an arena owned by the caller, which goes as a
// whole, so there is no partial state to unwind.
static bool try_read_cache(xtc_cache_t* cache, str_t path, md_allocator_i* alloc) {
    md_file_t file = {0};
    if (!md_run_cache_open(&file, &cache->header, path, MD_XTC_CACHE_MAGIC, MD_XTC_CACHE_VERSION)) {
        return false;
    }
    const size_t n = cache->header.num_frames;
    const bool ok =
        cache_read_block(file, (void**)&cache->frame_offsets, (n + 1) * sizeof(int64_t), alloc) &&
        cache_read_block(file, (void**)&cache->frame_times,   n * sizeof(double),        alloc) &&
        cache_read_block(file, (void**)&cache->frame_steps,   n * sizeof(int64_t),       alloc) &&
        cache_read_block(file, (void**)&cache->frame_boxes,   n * 9 * sizeof(float),     alloc) &&
        md_file_tell(file) == (int64_t)md_file_size(file);
    if (!ok) {
        MD_LOG_ERROR("XTC: the cache of '" STR_FMT "' is incomplete", STR_ARG(path));
    }
    md_file_close(&file);
    return ok;
}

static bool write_cache(const xtc_cache_t* cache, str_t path, const md_file_info_t* scanned) {
    md_file_t file = {0};
    const size_t n = cache->header.num_frames;
    if (!md_run_cache_create(&file, path, scanned, MD_XTC_CACHE_MAGIC, MD_XTC_CACHE_VERSION, cache->header.num_atoms, n)) {
        return false;
    }
    const bool ok =
        md_file_write(file, cache->frame_offsets, (n + 1) * sizeof(int64_t)) == (n + 1) * sizeof(int64_t) &&
        md_file_write(file, cache->frame_times,   n * sizeof(double))        == n * sizeof(double) &&
        md_file_write(file, cache->frame_steps,   n * sizeof(int64_t))       == n * sizeof(int64_t) &&
        md_file_write(file, cache->frame_boxes,   n * 9 * sizeof(float))     == n * 9 * sizeof(float);
    if (!ok) {
        MD_LOG_ERROR("XTC: failed to write the cache of '" STR_FMT "'", STR_ARG(path));
    }
    md_file_close(&file);
    return ok;
}

// What the file is, from its cache when that is current and from a scan otherwise (writing the
// cache unless told not to). path must be canonical; every array is allocated from alloc.
static bool xtc_index_load(xtc_cache_t* cache, str_t path, uint32_t flags, md_allocator_i* alloc) {
    MEMSET(cache, 0, sizeof(*cache));

    md_file_t file = {0};
    if (!md_file_open(&file, path, MD_FILE_READ)) {
        MD_LOG_ERROR("XTC: Failed to open file '" STR_FMT "'", STR_ARG(path));
        return false;
    }

    bool result = false;
    md_file_info_t scanned = {0};
    md_file_info_extract(file, &scanned);

    uint8_t frame_header_data[XTC_SMALL_HEADER_SIZE];
    md_xtc_header_t xtc_header = { 0 };
    if (md_file_read(file, frame_header_data, XTC_SMALL_HEADER_SIZE) != XTC_SMALL_HEADER_SIZE || !decode_header(frame_header_data, &xtc_header)) {
        MD_LOG_ERROR("XTC: Failed to read header of first frame, file may be corrupt or not a valid xtc trajectory");
        goto done;
    }
    if (xtc_header.natoms <= 0) {
        MD_LOG_ERROR("XTC: Number of atoms in trajectory was zero");
        goto done;
    }

    if (try_read_cache(cache, path, alloc)) {
        result = true;
        goto done;
    }

    MEMSET(cache, 0, sizeof(*cache));
    cache->header.num_atoms = xtc_header.natoms;

    md_array(int64_t) offsets = 0;
    md_array(double)  times   = 0;
    md_array(int64_t) steps   = 0;
    md_array(float)   boxes   = 0;
    cache->header.num_frames = xtc_scan(file, &offsets, &times, &steps, &boxes, alloc);
    if (!cache->header.num_frames || !offsets || !times) {
        MD_LOG_DEBUG("XTC: frame offsets or frame times was empty");
        goto done;
    }
    cache->frame_offsets = offsets;
    cache->frame_times   = times;
    cache->frame_steps   = steps;
    cache->frame_boxes   = boxes;

    if (!(flags & MD_RUN_FLAG_DISABLE_CACHE_WRITE)) {
        // If we fail to write the cache, that's ok, we can inform about it, but do not halt
        if (write_cache(cache, path, &scanned)) {
            MD_LOG_INFO("XTC: Successfully created cache file for '" STR_FMT "'", STR_ARG(path));
        }
    }
    result = true;

done:
    md_file_close(&file);
    return result;
}

// ### RUN ###

// <run>/atom/position. Everything it needs beyond the slice is in the table: the run is the path
// minus "/atom/position", and the file, offset and size of the frame are its source attributes.
static size_t xtc_position_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    const md_system_t* sys = (const md_system_t*)user_data;
    ASSERT(sys);
    const md_attributes_t* attributes = &sys->attributes;

    // A whole temporal virtual attribute is refused before it gets here, so the frame is fixed.
    if (!slice || slice->num_idx == 0 || slice->num_idx > 2) {
        return 0;
    }

    md_run_source_t src;
    if (!md_run_source(&src, attributes, attr, STR_LIT("atom/position"))) {
        return 0;
    }
    const int64_t* offsets = src.offset;
    const int64_t* sizes   = src.size;

    const uint32_t frame = slice->idx[0];
    const size_t num_atoms = attr->format.shape[1];
    const size_t frame_size = (size_t)sizes[frame];

    md_temp_scope_t temp = md_temp_begin();
    size_t written = 0;

    // The decoder reads the bit stream a word at a time and may look past the end of the frame.
    uint8_t* frame_data = md_temp_alloc(temp, ALIGN_TO(frame_size, 16) + MD_XTC_STREAM_GUARD_BYTES);
    if (frame_data) {
        MEMSET(frame_data + frame_size, 0, ALIGN_TO(frame_size, 16) + MD_XTC_STREAM_GUARD_BYTES - frame_size);
    }
    // Through io: inside an extraction context the file stays open from frame to frame, which is
    // what makes streaming from a cluster's file server bearable. Without one it is opened here.
    const str_t file_path = src.path;
    const size_t read = frame_data ? md_attribute_io_read_at(io, file_path, offsets[frame], frame_data, frame_size) : 0;
    if (frame_data && read != frame_size) {
        MD_LOG_ERROR("XTC: Failed to read frame %u from '" STR_FMT "', expected %zu bytes, got %zu", frame, STR_ARG(file_path), frame_size, read);
    }
    if (frame_data && read == frame_size) {
        if (slice->num_idx == 1 && cap == num_atoms * 3) {
            // The whole frame, which is the case that matters: decoded straight into the caller.
            if (xtc_decode_frame_data_scaled(frame_data, frame_size, NULL, (float*)dst, num_atoms, 10.0f)) {
                written = cap;
            }
        } else if (slice->num_idx == 2 && cap == 3 && slice->idx[1] < num_atoms) {
            // One atom: the stream is not seekable, so the frame is decoded and one value kept.
            float* xyz = md_temp_alloc(temp, num_atoms * 3 * sizeof(float));
            if (xyz && xtc_decode_frame_data_scaled(frame_data, frame_size, NULL, xyz, num_atoms, 10.0f)) {
                MEMCPY(dst, xyz + (size_t)slice->idx[1] * 3, 3 * sizeof(float));
                written = cap;
            }
        }
    }

    md_temp_end(temp);
    return written;
}

bool md_xtc_system_publish_run(md_system_t* sys, str_t filename, str_t run, uint32_t flags) {
    ASSERT(sys);
    if (str_empty(run)) {
        MD_LOG_ERROR("XTC: no run to publish into");
        return false;
    }
    char path_buf[4096];
    const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
    const str_t path = {path_buf, path_len};

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    bool result = false;

    xtc_cache_t index;
    if (!xtc_index_load(&index, path, flags, arena)) {
        goto done;
    }

    const size_t F = index.header.num_frames;
    const size_t N = index.header.num_atoms;

    int64_t* sizes = md_alloc(arena, F * sizeof(int64_t));
    float*   boxes = md_alloc(arena, F * 9 * sizeof(float));
    for (size_t i = 0; i < F; ++i) {
        sizes[i] = index.frame_offsets[i + 1] - index.frame_offsets[i];
    }
    for (size_t i = 0; i < F * 9; ++i) {
        boxes[i] = index.frame_boxes[i] * 10.0f;   // nm to Angstrom, as the coordinates
    }

    const md_attribute_virtual_t virt = { .provider = xtc_position_provider, .user_data = sys };
    const md_run_desc_t desc = {
        .num_frames    = F,
        .num_atoms     = N,
        .time          = index.frame_times,
        .time_unit     = md_unit_picosecond(),
        .step          = index.frame_steps,
        .unitcell      = boxes,
        .source_path   = path,
        .source_offset = index.frame_offsets,
        .source_size   = sizes,
        .position_virt = &virt,
    };
    result = md_run_publish(sys, run, &desc);

done:
    md_arena_allocator_destroy(arena);
    return result;
}

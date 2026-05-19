#include <md_xtc.h>

#include <md_system.h>
#include <md_trajectory.h>

#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_vec_math.h>

#if defined(__AVX2__)
#ifndef LIBDIVIDE_AVX2
#define LIBDIVIDE_AVX2
#endif
#define MD_XTC_AVX2_FASTPATH 1
#else
#define MD_XTC_AVX2_FASTPATH 0
#endif

#include <libdivide.h>
#include <stdio.h>

#if defined(__SIZEOF_INT128__)
#define HAS_INT128_T
#endif

#define MD_XTC_CACHE_MAGIC   0x8281237612371
#define MD_XTC_CACHE_VERSION 3
#define MD_XTC_TRAJ_MAGIC 0x162365dac721995
#define MD_XTC_TRAJ_READER_MAGIC 0x162365dac721996

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

#if MD_XTC_AVX2_FASTPATH
#define DIV32_T struct libdivide_u32_branchfree_t
#define DIV32_INIT(x) libdivide_u32_branchfree_gen(x)
#define DIV32_VEC(x, y) libdivide_u32_branchfree_do_vec256(x, y)
#endif

typedef md_128i v4i_t;
#define v4i_set(x, y, z, w) md_mm_set_epi32(w, z, y, x)
#define v4i_set1(x)         md_mm_set1_epi32(x)
#define v4i_add(a, b)       md_mm_add_epi32(a, b)
#define v4i_sub(a, b)       md_mm_sub_epi32(a, b)
#define v4i_load(addr)      md_mm_loadu_epi32(addr)

// Trajectory opaque data for xtc
typedef struct xtc_t {
    uint64_t magic;
	str_t filepath;
    md_array(int64_t) frame_offsets;
    md_trajectory_header_t header;
    md_allocator_i* alloc;
} xtc_t;

// Trajectory reader state
typedef struct xtc_reader_t {
	uint64_t magic;
    md_file_t file;
    md_array(uint8_t) frame_data; // scratch buffer for compressed frame data, resized as needed for each frame
    size_t num_atoms;
    size_t num_frames;
	const int64_t* frame_offsets; // pointer into xtc_t.frame_offsets for quick access
    md_allocator_i*   arena;
} xtc_reader_t;

// Number of guard bytes that must be readable past the end of a bitstream so
// the stateless unaligned reads in extract_bits_be_raw_* never overrun.
#define MD_XTC_STREAM_GUARD_BYTES 64
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

#if MD_XTC_AVX2_FASTPATH
typedef struct unpack_data32_t {
    uint32_t   size_y;
    uint32_t   size_z;
    bit_data_t bit;
    DIV32_T    div_zy;
    DIV32_T    div_z;
} unpack_data32_t;
#endif

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
// All scalar extractors load up to 16 bytes from `base + (bit_offset >> 3)`;
// the AVX2 small-run extractor reads two 32-byte chunks. The
// caller MUST ensure that `MD_XTC_STREAM_GUARD_BYTES` bytes past the
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

#if MD_XTC_AVX2_FASTPATH
static FORCE_INLINE __m256i bswap_epi32x8_avx2(__m256i v) {
    const __m256i shuf = _mm256_setr_epi8(
        3,  2,  1,  0,  7,  6,  5,  4, 11, 10,  9,  8, 15, 14, 13, 12,
        3,  2,  1,  0,  7,  6,  5,  4, 11, 10,  9,  8, 15, 14, 13, 12);
    return _mm256_shuffle_epi8(v, shuf);
}

static FORCE_INLINE unpack_data32_t init_small_unpack32_avx2(uint32_t smallidx) {
    ASSERT(FIRSTIDX <= smallidx && smallidx <= 32);
    const uint32_t smallsize = magicints[smallidx];
    unpack_data32_t unpack = {
        .size_y = smallsize,
        .size_z = smallsize,
        .div_zy = DIV32_INIT(smallsize * smallsize),
        .div_z  = DIV32_INIT(smallsize),
    };
    init_unpack_bit_data(&unpack.bit, smallidx);
    return unpack;
}

static FORCE_INLINE __m256i extract_bits_be_raw_32x8_avx2(const uint8_t* base, size_t bit_offset, uint32_t bit_length) {
    ASSERT(0 < bit_length && bit_length <= 32);

    const uint8_t* ptr = base + (bit_offset >> 3);
    const uint32_t bit_in_byte = (uint32_t)(bit_offset & 7);

    __m256i w0 = _mm256_loadu_si256((const __m256i*)((const void*)ptr));
    __m256i w1 = _mm256_loadu_si256((const __m256i*)((const void*)(ptr + 32)));
    w0 = bswap_epi32x8_avx2(w0);
    w1 = bswap_epi32x8_avx2(w1);

    const __m256i lane = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    const __m256i start = _mm256_add_epi32(_mm256_mullo_epi32(lane, _mm256_set1_epi32((int)bit_length)), _mm256_set1_epi32((int)bit_in_byte));
    const __m256i word_idx = _mm256_srli_epi32(start, 5);
    const __m256i bit_idx  = _mm256_and_si256(start, _mm256_set1_epi32(31));

    const __m256i word = _mm256_permutevar8x32_epi32(w0, word_idx);
    const __m256i next_idx = _mm256_add_epi32(word_idx, _mm256_set1_epi32(1));
    const __m256i next0 = _mm256_permutevar8x32_epi32(w0, next_idx);
    const __m256i next1 = _mm256_permutevar8x32_epi32(w1, _mm256_sub_epi32(next_idx, _mm256_set1_epi32(8)));
    const __m256i next = _mm256_blendv_epi8(next0, next1, _mm256_cmpgt_epi32(next_idx, _mm256_set1_epi32(7)));

    const __m256i lo = _mm256_sllv_epi32(word, bit_idx);
    const __m256i hi = _mm256_srlv_epi32(next, _mm256_sub_epi32(_mm256_set1_epi32(32), bit_idx));
    return _mm256_srlv_epi32(_mm256_or_si256(lo, hi), _mm256_set1_epi32((int)(32 - bit_length)));
}

static FORCE_INLINE void unpack_coord32x8_avx2(__m256i* out_x, __m256i* out_y, __m256i* out_z, __m256i w, const unpack_data32_t* unpack) {
    const __m256i t0 = _mm256_sllv_epi32(w, _mm256_set1_epi32((int)unpack->bit.sml_shift));
    const __m256i t1 = _mm256_and_si256(t0, _mm256_set1_epi32((int)~0xFFu));
    const __m256i t2 = _mm256_and_si256(w,  _mm256_set1_epi32((int)unpack->bit.part_mask));
    const __m256i t  = bswap_epi32x8_avx2(_mm256_or_si256(t1, t2));
    const uint32_t packed_bits = (unpack->bit.num_of_bits + 7u) & ~7u;
    const __m256i v = _mm256_srlv_epi32(t, _mm256_set1_epi32((int)(32u - packed_bits)));

    const __m256i x  = DIV32_VEC(v,  &unpack->div_zy);
    const __m256i yz = DIV32_VEC(v,  &unpack->div_z);
    const __m256i y  = _mm256_sub_epi32(yz, _mm256_mullo_epi32(x,  _mm256_set1_epi32((int)unpack->size_y)));
    const __m256i z  = _mm256_sub_epi32(v,  _mm256_mullo_epi32(yz, _mm256_set1_epi32((int)unpack->size_z)));

    *out_x = x;
    *out_y = y;
    *out_z = z;
}

static FORCE_INLINE __m256i prefix_sum_epi32x8_avx2(__m256i v) {
    v = _mm256_add_epi32(v, _mm256_slli_si256(v, 4));
    v = _mm256_add_epi32(v, _mm256_slli_si256(v, 8));
    const __m256i lo_sum = _mm256_permutevar8x32_epi32(v, _mm256_set1_epi32(3));
    return _mm256_add_epi32(v, _mm256_blend_epi32(_mm256_setzero_si256(), lo_sum, 0xF0));
}

static FORCE_INLINE __m256i propagate_small_axis8_avx2(__m256i coord, int32_t prev, int32_t smallnum) {
    const __m256i delta = _mm256_sub_epi32(coord, _mm256_set1_epi32(smallnum));
    return _mm256_add_epi32(prefix_sum_epi32x8_avx2(delta), _mm256_set1_epi32(prev));
}

static FORCE_INLINE void decode_small_run8_avx2(__m256i* out_x, __m256i* out_y, __m256i* out_z, const uint8_t* stream, size_t bit_offset, const unpack_data32_t* unpack, int32_t prev_x, int32_t prev_y, int32_t prev_z, int32_t smallnum) {
    __m256i x, y, z;
    __m256i w = extract_bits_be_raw_32x8_avx2(stream, bit_offset, unpack->bit.num_of_bits);
    unpack_coord32x8_avx2(&x, &y, &z, w, unpack);
    *out_x = propagate_small_axis8_avx2(x, prev_x, smallnum);
    *out_y = propagate_small_axis8_avx2(y, prev_y, smallnum);
    *out_z = propagate_small_axis8_avx2(z, prev_z, smallnum);
}

static FORCE_INLINE void store_small_run8_axis_soa_avx2(float* RESTRICT out, int atom_idx, int run_count, __m256i coord, int32_t big, __m256 scale_vec, float scale) {
    out[atom_idx] = (float)_mm_cvtsi128_si32(_mm256_castsi256_si128(coord)) * scale;

    __m256i ordered = _mm256_blend_epi32(coord, _mm256_set1_epi32(big), 0x01);
    __m256i mask = _mm256_cmpgt_epi32(_mm256_set1_epi32(run_count), _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7));
    _mm256_maskstore_ps(out + atom_idx + 1, mask, _mm256_mul_ps(_mm256_cvtepi32_ps(ordered), scale_vec));
}

static FORCE_INLINE void store_small_run8_soa_avx2(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z, int atom_idx, int run_count, __m256i x, __m256i y, __m256i z, int32_t big_x, int32_t big_y, int32_t big_z, __m256 scale_vec, float scale) {
    store_small_run8_axis_soa_avx2(out_x, atom_idx, run_count, x, big_x, scale_vec, scale);
    store_small_run8_axis_soa_avx2(out_y, atom_idx, run_count, y, big_y, scale_vec, scale);
    store_small_run8_axis_soa_avx2(out_z, atom_idx, run_count, z, big_z, scale_vec, scale);
}

static FORCE_INLINE int32_t extract_lane_epi32x8_avx2(__m256i v, int lane) {
    __m256i splat = _mm256_permutevar8x32_epi32(v, _mm256_set1_epi32(lane));
    return _mm_cvtsi128_si32(_mm256_castsi256_si128(splat));
}

static FORCE_INLINE v4i_t last_coord_avx2(__m256i x, __m256i y, __m256i z, int run_count) {
    int lane = run_count - 1;
    return v4i_set(extract_lane_epi32x8_avx2(x, lane), extract_lane_epi32x8_avx2(y, lane), extract_lane_epi32x8_avx2(z, lane), 0);
}
#endif

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

static inline void extract_int32(int32_t* out, size_t count, const uint8_t* ptr) {
	MEMCPY(out, ptr, sizeof(int32_t) * count);
#if __LITTLE_ENDIAN__
	for (size_t i = 0; i < count; ++i) {
		out[i] = BSWAP32(out[i]);
	}
#endif
}

static inline void extract_float(float* out, size_t count, const uint8_t* ptr) {
    for (size_t i = 0; i < count; ++i) {
        int32_t ival;
        MEMCPY(&ival, ptr + i * sizeof(int32_t), sizeof(int32_t));
#if __LITTLE_ENDIAN__
        ival = BSWAP32(ival);
#endif
        float fval;
        MEMCPY(&fval, &ival, sizeof(float));
        out[i] = fval;
    }
}

static inline bool decode_header(const uint8_t* frame_ptr, md_xtc_header_t* out_header) {
    // Extract header
    int magic;
    extract_int32(&magic,  1, frame_ptr);
    if (magic != XTC_MAGIC) {
        MD_LOG_ERROR("XTC: Magic number did not match");
        return false;
	}

    extract_int32(&out_header->natoms, 1, frame_ptr + 4);
    extract_int32(&out_header->step,   1, frame_ptr + 8);
	extract_float(&out_header->time,   1, frame_ptr + 12);
	extract_float((float*)out_header->box, 9, frame_ptr + 16);

    return true;
}

size_t md_xtc_read_frame_offsets_and_times(md_file_t xdr, md_array(int64_t)* frame_offsets, md_array(double)* frame_times, md_allocator_i* alloc) {
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

        if (!md_file_read(xdr, &framebytes, sizeof(int32_t))) {
            MD_LOG_ERROR("XTC: Failed to read framebytes");
            return 0;
        }
#if __LITTLE_ENDIAN__
        framebytes = BSWAP32(framebytes);
#endif
        framebytes = (framebytes + 3) & ~0x03; /* Rounding to the next 32-bit boundary */

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
            if (!md_file_read(xdr, &framebytes, sizeof(int32_t))) {
                MD_LOG_ERROR("XTC: Failed to read framebytes");
                goto done;
            }
#if __LITTLE_ENDIAN__
            framebytes = BSWAP32(framebytes);
#endif
            framebytes = (framebytes + 3) & ~0x03; /* Rounding to the next 32-bit boundary */

            /* Skip `framebytes` to next header */
            if (!md_file_seek(xdr, framebytes, MD_FILE_CUR)) {
                MD_LOG_DEBUG("XTC: encountered corrupted frame");
                goto done;
            }

            /* Store position in `offsets`, adjust for header */
            md_array_push(*frame_offsets, offset, alloc);
            md_array_push(*frame_times, xtc_header.time, alloc);
            num_frames += 1;
        }
        // Add last offset
        md_array_push(*frame_offsets, filesize, alloc);
    }
done:
    return num_frames;
}

bool md_xtc_decode_frame_data(const uint8_t* frame_ptr, size_t frame_bytes, md_xtc_header_t* out_header, float* out_coords, size_t num_atoms) {
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
	}

    if (!out_coords) {
        return true;
    }

    int natoms;
	extract_int32(&natoms, 1, frame_ptr + XTC_SMALL_HEADER_SIZE - 4);

    if (natoms != (int)num_atoms) {
        MD_LOG_ERROR("XTC: Number of atoms in frame header does not match expected number of atoms");
        return false;
    }

    size_t offset = XTC_SMALL_HEADER_SIZE;
    if (natoms <= 9) {
		// No compression for 9 atoms or less, just read the coordinates directly
		extract_float(out_coords, (size_t)natoms * 3, frame_ptr + offset);
        return true;
    }

    float precision;
    int32_t minint[3], maxint[3], smallidx;
	extract_float(&precision, 1, frame_ptr + offset); offset += 4;
	extract_int32(minint,     3, frame_ptr + offset); offset += 12;
	extract_int32(maxint,     3, frame_ptr + offset); offset += 12;
	extract_int32(&smallidx,  1, frame_ptr + offset); offset += 4;

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
    extract_int32(&num_bytes, 1, frame_ptr + offset); offset += 4;
    (void)num_bytes;

    const uint8_t* stream = frame_ptr + offset;
    size_t bit_offset = 0;

    float* lfp = out_coords;
    md_128 invp = md_mm_set1_ps(1.0f / precision);
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
    extract_int32(&natoms, 1, frame_ptr + XTC_SMALL_HEADER_SIZE - 4);

    if (natoms != (int)num_atoms) {
        MD_LOG_ERROR("XTC: Number of atoms in frame header does not match expected number of atoms");
        return false;
    }

    size_t offset = XTC_SMALL_HEADER_SIZE;
    if (natoms <= 9) {
        for (int i = 0; i < natoms; ++i) {
            float coord[3];
            extract_float(coord, 3, frame_ptr + offset + (size_t)i * XTC_SMALL_COORDS_SIZE);
            out_x[i] = coord[0] * scale;
            out_y[i] = coord[1] * scale;
            out_z[i] = coord[2] * scale;
        }
        return true;
    }

    float precision;
    int32_t minint[3], maxint[3], smallidx;
    extract_float(&precision, 1, frame_ptr + offset); offset += 4;
    extract_int32(minint,     3, frame_ptr + offset); offset += 12;
    extract_int32(maxint,     3, frame_ptr + offset); offset += 12;
    extract_int32(&smallidx,  1, frame_ptr + offset); offset += 4;

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

#if MD_XTC_AVX2_FASTPATH
    unpack_data32_t sml_unpack32 = {0};
    if (smallidx <= 32) {
        sml_unpack32 = init_small_unpack32_avx2((uint32_t)smallidx);
    }
#endif

    int32_t num_bytes = 0;
    extract_int32(&num_bytes, 1, frame_ptr + offset); offset += 4;
    (void)num_bytes;

    const uint8_t* stream = frame_ptr + offset;
    size_t bit_offset = 0;

    const float coord_scale_val = scale / precision;
    md_128 coord_scale = md_mm_set1_ps(coord_scale_val);
#if MD_XTC_AVX2_FASTPATH
    __m256 coord_scale_avx = _mm256_set1_ps(coord_scale_val);
#endif
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

#if MD_XTC_AVX2_FASTPATH
            if (sml_bits <= 32) {
                ASSERT(run_count <= 8); // 8 runs of 32 bits fits in 256-bit register
                const int32_t big_x = simde_mm_extract_epi32(prevcoord, 0);
                const int32_t big_y = simde_mm_extract_epi32(prevcoord, 1);
                const int32_t big_z = simde_mm_extract_epi32(prevcoord, 2);

                __m256i x, y, z;
                decode_small_run8_avx2(&x, &y, &z, stream, bit_offset, &sml_unpack32, big_x, big_y, big_z, smallnum);
                bit_offset += (size_t)sml_bits * run_count;

                store_small_run8_soa_avx2(out_x, out_y, out_z, atom_idx, run_count, x, y, z, big_x, big_y, big_z, coord_scale_avx, coord_scale_val);
                atom_idx += run_count + 1;
                thiscoord = last_coord_avx2(x, y, z, run_count);
            } else
#endif
            {
                uint64_t w = extract_bits_be_raw_57(stream, bit_offset, sml_bits);
                v4i_t coord = unpack_coord64(w, &sml_unpack);
                bit_offset += sml_bits;
                thiscoord = v4i_add(coord, v4i_sub(thiscoord, vsmall));

                // Write second before first coord as is required by the compression scheme
                write_coord_soa(out_x, out_y, out_z, atom_idx++, thiscoord, coord_scale);
                write_coord_soa(out_x, out_y, out_z, atom_idx++, prevcoord, coord_scale);

                for (int i = 1; i < run_count; ++i) {
                    w = extract_bits_be_raw_57(stream, bit_offset, sml_bits);
                    coord = unpack_coord64(w, &sml_unpack);
                    bit_offset += sml_bits;
                    thiscoord = v4i_add(coord, v4i_sub(thiscoord, vsmall));
                    write_coord_soa(out_x, out_y, out_z, atom_idx++, thiscoord, coord_scale);
                }
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
#if MD_XTC_AVX2_FASTPATH
            if (smallidx <= 32) {
                sml_unpack32 = init_small_unpack32_avx2((uint32_t)smallidx);
            }
#endif
        }
    }

done:
    return atom_idx == natoms;
}

bool md_xtc_decode_frame_data_soa(const uint8_t* frame_ptr, size_t frame_bytes, md_xtc_header_t* out_header, float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z, size_t num_atoms) {
    return md_xtc_decode_frame_data_soa_scaled(frame_ptr, frame_bytes, out_header, out_x, out_y, out_z, num_atoms, 1.0f);
}

static bool xtc_get_header(struct md_trajectory_o* inst, md_trajectory_header_t* header) {
    xtc_t* xtc = (xtc_t*)inst;
    ASSERT(xtc);
    ASSERT(xtc->magic == MD_XTC_TRAJ_MAGIC);
    ASSERT(header);

    *header = xtc->header;
    return true;
}

static bool xtc_reader_load_frame(struct md_trajectory_reader_o* inst, int64_t frame_idx, md_trajectory_frame_header_t* out_header, float* out_x, float* out_y, float* out_z) {
    ASSERT(inst);

    xtc_reader_t* xtc = (xtc_reader_t*)inst;
    bool result = false;
    if (md_file_valid(xtc->file)) {
        if (!xtc->frame_offsets) {
            MD_LOG_ERROR("XTC: Frame offsets is empty");
            return 0;
        }

        if (frame_idx < 0 || (int64_t)xtc->num_frames <= frame_idx) {
            MD_LOG_ERROR("XTC: Frame index is out of range");
            return 0;
        }
		const int64_t beg = xtc->frame_offsets[frame_idx];
		const int64_t end = xtc->frame_offsets[frame_idx + 1];
        if (end <= beg) {
            MD_LOG_ERROR("XTC: Invalid frame offset range");
            return false;
		}
		const size_t frame_size = end - beg;

        md_array_ensure(xtc->frame_data, ALIGN_TO(frame_size, 16) + MD_XTC_STREAM_GUARD_BYTES, xtc->arena);
        size_t read_size = md_file_read_at(xtc->file, beg, xtc->frame_data, frame_size);

        if (read_size == frame_size) {
            md_xtc_header_t xtc_header = {0};

            if (md_xtc_decode_frame_data_soa_scaled(xtc->frame_data, frame_size, &xtc_header, out_x, out_y, out_z, xtc->num_atoms, 10.0f)) {
                if (out_header) {
                    out_header->num_atoms = xtc_header.natoms;
                    out_header->index     = xtc_header.step;
                    out_header->timestamp = xtc_header.time;
                    out_header->unitcell  = md_unitcell_from_matrix_float(xtc_header.box);
                }

                result = true;
            } else {
                MD_LOG_ERROR("XTC: Failed to decode frame data");
            }
        } else {
            MD_LOG_ERROR("XTC: Failed to read frame data from file, expected %zu bytes, got %zu bytes", frame_size, read_size);
        }
    }

    return result;
}

static void xtc_trajectory_reader_free(struct md_trajectory_reader_i* reader) {
    if (!reader) {
        return;
    }
    xtc_reader_t* inst = (xtc_reader_t*)reader->inst;
    if (inst) {
        ASSERT(inst->magic == MD_XTC_TRAJ_READER_MAGIC);
        if (md_file_valid(inst->file)) {
            md_file_close(&inst->file);
        }
        md_arena_allocator_destroy(inst->arena);
    }
    MEMSET(reader, 0, sizeof(md_trajectory_reader_i));
}

static bool xtc_trajectory_reader_init(md_trajectory_reader_i* reader, struct md_trajectory_o* traj_inst) {
    ASSERT(reader);
    ASSERT(traj_inst);

    xtc_t* xtc = (xtc_t*)traj_inst;
    ASSERT(xtc->magic == MD_XTC_TRAJ_MAGIC);

    md_file_t file = {0};
    if (!md_file_open(&file, xtc->filepath, MD_FILE_READ)) {
        MD_LOG_ERROR("XTC: Failed to open file '" STR_FMT "'", STR_ARG(xtc->filepath));
        return false;
    }

    MEMSET(reader, 0, sizeof(md_trajectory_reader_i));

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    xtc_reader_t* inst = md_alloc(arena, sizeof(xtc_reader_t));
    MEMSET(inst, 0, sizeof(xtc_reader_t));
    inst->magic = MD_XTC_TRAJ_READER_MAGIC;
    inst->file = file;
    inst->arena = arena;
    inst->num_atoms = xtc->header.num_atoms;
    inst->num_frames = xtc->header.num_frames;
    inst->frame_offsets = xtc->frame_offsets;

    reader->inst = (struct md_trajectory_reader_o*)inst;
    reader->load_frame = xtc_reader_load_frame;
    reader->free = xtc_trajectory_reader_free;

    return true;
}

static bool xtc_load_frame(struct md_trajectory_o* inst, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(inst);

    xtc_t* xtc = (xtc_t*)inst;
    if (xtc->magic != MD_XTC_TRAJ_MAGIC) {
        MD_LOG_ERROR("XTC: Error when decoding frame coord, xtc magic did not match");
        return false;
    }
    if (!xtc->frame_offsets) {
        MD_LOG_ERROR("XTC: Frame offsets is empty");
        return 0;
    }

    if (frame_idx < 0 || (int64_t)xtc->header.num_frames <= frame_idx) {
        MD_LOG_ERROR("XTC: Frame index is out of range");
        return 0;
    }

    bool result = false;
	md_file_t file = { 0 };
    if (md_file_open(&file, xtc->filepath, MD_FILE_READ)) {

        const int64_t beg = xtc->frame_offsets[frame_idx];
        const int64_t end = xtc->frame_offsets[frame_idx + 1];
        if (end <= beg) {
            MD_LOG_ERROR("XTC: Invalid frame offset range");
            return false;
        }
        const size_t frame_size = end - beg;
		const size_t alloc_size = ALIGN_TO(frame_size, 16) + MD_XTC_STREAM_GUARD_BYTES;

		uint8_t* frame_data = calloc(alloc_size, 1);
        if (!frame_data) {
            MD_LOG_ERROR("XTC: Failed to allocate memory for frame data");
            return false;
		}

        size_t read_size = md_file_read_at(file, beg, frame_data, frame_size);
		if (read_size == frame_size) {
            size_t num_atoms = xtc->header.num_atoms;
            md_xtc_header_t xtc_header = {0};

            if (md_xtc_decode_frame_data_soa_scaled(frame_data, frame_size, &xtc_header, x, y, z, num_atoms, 10.0f)) {
                if (header) {
                    for (int i = 0; i < 3; ++i) {
                        xtc_header.box[i][0] *= 10.0f;
                        xtc_header.box[i][1] *= 10.0f;
                        xtc_header.box[i][2] *= 10.0f;
			        }
                    header->num_atoms = xtc_header.natoms;
			        header->index     = xtc_header.step;
			        header->timestamp = xtc_header.time;
			        header->unitcell  = md_unitcell_from_matrix_float(xtc_header.box);
                }

                result = true;
            } else {
                MD_LOG_ERROR("XTC: Failed to decode frame data");
            }
        } else {
			MD_LOG_ERROR("XTC: Failed to read frame data from file, expected %zu bytes, got %zu bytes", frame_size, read_size);
        }

		free(frame_data);
		md_file_close(&file);
    }

    return result;
}

typedef struct xtc_cache_t {
    md_trajectory_cache_header_t header;
    int64_t* frame_offsets;
    double*  frame_times;
} xtc_cache_t;

static bool try_read_cache(xtc_cache_t* cache, str_t cache_file, size_t traj_num_bytes, md_allocator_i* alloc) {
    ASSERT(cache);
    ASSERT(alloc);

    bool result = false;
    md_file_t file = {0};
    if (md_file_open(&file, cache_file, MD_FILE_READ)) {
        if (md_file_read(file, &cache->header, sizeof(cache->header)) != sizeof(cache->header)) {
            MD_LOG_ERROR("XTC trajectory cache: failed to read header");
            goto done;
        }

        if (cache->header.magic != MD_XTC_CACHE_MAGIC) {
            MD_LOG_ERROR("XTC trajectory cache: magic was incorrect or corrupt");
            goto done;
        }
        if (cache->header.version != MD_XTC_CACHE_VERSION) {
            MD_LOG_INFO("XTC trajectory cache: version mismatch, expected %i, got %i", MD_XTC_CACHE_VERSION, (int)cache->header.version);
        }
        if (cache->header.num_bytes != traj_num_bytes) {
            MD_LOG_INFO("XTC trajectory cache: trajectory size mismatch, expected %zu, got %zu", traj_num_bytes, cache->header.num_bytes);
        }
        if (cache->header.num_atoms == 0) {
            MD_LOG_ERROR("XTC trajectory cache: num atoms was zero");
            goto done;
        }
        if (cache->header.num_frames == 0) {
            MD_LOG_ERROR("XTC trajectory cache: num frames was zero");
            goto done;
        }

        const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
        cache->frame_offsets = md_alloc(alloc, offset_bytes);
        if (md_file_read(file, cache->frame_offsets, offset_bytes) != offset_bytes) {
            MD_LOG_ERROR("XTC trajectory cache: Failed to read offset data");
            md_free(alloc, cache->frame_offsets, offset_bytes);
            goto done;
        }

        const size_t time_bytes = cache->header.num_frames * sizeof(double);
        cache->frame_times = md_alloc(alloc, time_bytes);
        if (md_file_read(file, cache->frame_times, time_bytes) != time_bytes) {
        	MD_LOG_ERROR("XTC trajectory cache: times are incomplete");
        	md_free(alloc, cache->frame_offsets, offset_bytes);
        	md_free(alloc, cache->frame_times, time_bytes);
        	goto done;
        }

        // Test position in file, we expect to be at the end of the file
        if (md_file_tell(file) != (int64_t)md_file_size(file)) {
        	MD_LOG_ERROR("XTC trajectory cache: file position was not at the end of the file");
        	md_free(alloc, cache->frame_offsets, offset_bytes);
        	md_free(alloc, cache->frame_times, time_bytes);
        	goto done;
        }

        result = true;
    done:
        md_file_close(&file);
    }
    return result;
}

static bool write_cache(const xtc_cache_t* cache, str_t cache_file) {
    bool result = false;

    md_file_t file = {0};
    if (!md_file_open(&file, cache_file, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        MD_LOG_INFO("XTC trajectory cache: could not open file '"STR_FMT"'", STR_ARG(cache_file));
        return false;
    }

    if (md_file_write(file, &cache->header, sizeof(cache->header)) != sizeof(cache->header)) {
        MD_LOG_ERROR("XTC trajectory cache: failed to write header");
        goto done;
    }

    const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
    if (md_file_write(file, cache->frame_offsets, offset_bytes) != offset_bytes) {
        MD_LOG_ERROR("Failed to write offset cache, offsets");
        goto done;
    }

    const size_t time_bytes = cache->header.num_frames * sizeof(double);
    if (md_file_write(file, cache->frame_times, time_bytes) != time_bytes) {
    	MD_LOG_ERROR("Failed to write offset cache, times");
    	goto done;
    }

    result = true;

done:
    md_file_close(&file);
    return result;
}

static void md_xtc_trajectory_free(md_trajectory_i* traj) {
    ASSERT(traj);
    ASSERT(traj->inst);
    xtc_t* xtc = (xtc_t*)traj->inst;
    if (xtc->magic != MD_XTC_TRAJ_MAGIC) {
        MD_LOG_ERROR("XTC: Cannot free trajectory, is not a valid XTC trajectory.");
        ASSERT(false);
        return;
    }
    MEMSET(traj, 0, sizeof(md_trajectory_i));
    md_arena_allocator_destroy(xtc->alloc);
}

md_trajectory_i* md_xtc_trajectory_create(str_t filename, md_allocator_i* ext_alloc, uint32_t flags) {
    ASSERT(ext_alloc);
    md_allocator_i* alloc = md_arena_allocator_create(ext_alloc, MEGABYTES(1));

    char path_buf[4096];
    size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
    str_t path = {path_buf, path_len};

    md_file_t file = {0};
    if (md_file_open(&file, path, MD_FILE_READ)) {
        const size_t filesize = md_file_size(file);

        uint8_t frame_header_data[XTC_SMALL_HEADER_SIZE];
        if (md_file_read(file, frame_header_data, XTC_SMALL_HEADER_SIZE) != XTC_SMALL_HEADER_SIZE) {
            MD_LOG_ERROR("XTC: Failed to read header of first frame, file may be corrupt or not a valid xtc trajectory");
            goto fail;
        }

        md_xtc_header_t xtc_header = { 0 };
        if (!decode_header(frame_header_data, &xtc_header)) {
            MD_LOG_ERROR("XTC: Failed to decode header of first frame, file may be corrupt or not a valid xtc trajectory");
            goto fail;
        }

        if (xtc_header.natoms == 0) {
            MD_LOG_ERROR("XTC: Number of atoms in trajectory was zero");
            goto fail;
        }

        char cache_buf[4096];
		int len = snprintf(cache_buf, sizeof(cache_buf), STR_FMT ".cache", STR_ARG(path));
		ASSERT(0 < len && len < (int)sizeof(cache_buf));
        str_t cache_path = { cache_buf, (size_t)len };

        xtc_cache_t cache = {0};
        if (!try_read_cache(&cache, cache_path, filesize, alloc)) {
            cache.header.magic = MD_XTC_CACHE_MAGIC;
            cache.header.version = MD_XTC_CACHE_VERSION;
            cache.header.num_bytes = filesize;
            cache.header.num_atoms = xtc_header.natoms;
            cache.header.num_frames = md_xtc_read_frame_offsets_and_times(file, &cache.frame_offsets, &cache.frame_times, alloc);
            if (!cache.header.num_frames) {
                goto fail;
            }
            if (!cache.frame_offsets || !cache.frame_times) {
                MD_LOG_DEBUG("XTC: frame offsets or frame times was empty");
                goto fail;
            }

            if (!(flags & MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE)) {
                // If we fail to write the cache, that's ok, we can inform about it, but do not halt
                if (write_cache(&cache, cache_path)) {
                    MD_LOG_INFO("XTC: Successfully created cache file for '" STR_FMT "'", STR_ARG(path));
                }
            }
        }

        void* mem = md_alloc(alloc, sizeof(md_trajectory_i) + sizeof(xtc_t));
        ASSERT(mem);
        MEMSET(mem, 0, sizeof(md_trajectory_i) + sizeof(xtc_t));

        md_trajectory_i* traj = mem;
        xtc_t* xtc = (xtc_t*)(traj + 1);

        xtc->magic = MD_XTC_TRAJ_MAGIC;
        xtc->filepath = str_copy(path, alloc);
        xtc->alloc = alloc;
        xtc->frame_offsets = cache.frame_offsets;

        xtc->header = (md_trajectory_header_t) {
            .num_frames = cache.header.num_frames,
            .num_atoms = xtc_header.natoms,
            .time_unit = md_unit_picosecond(),
            .frame_times = cache.frame_times,
        };

        traj->inst = (struct md_trajectory_o*)xtc;
        traj->free = md_xtc_trajectory_free;
        traj->get_header = xtc_get_header;
		traj->init_reader = xtc_trajectory_reader_init;

        md_file_close(&file);
        return traj;
    }
fail:
    if (md_file_valid(file)) md_file_close(&file);
    md_arena_allocator_destroy(alloc);
    return NULL;
}

// Attach convenience wrapper: create trajectory and attach to system
bool md_xtc_attach_from_file(struct md_system_t* sys, str_t filename, uint32_t flags) {
    if (!sys) return false;
    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }
    md_trajectory_i* traj = md_xtc_trajectory_create(filename, sys->alloc, flags);
    if (!traj) return false;
    md_system_attach_trajectory(sys, traj);
    return true;
}

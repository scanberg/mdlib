#include <md_xtc.h>

#include <md_util.h>
#include <md_trajectory.h>

#include <core/md_platform.h>
#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_str_builder.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_vec_math.h>

#include <libdivide.h>

#include <string.h>
#include <stdio.h>

#if defined(__SIZEOF_INT128__)
#define HAS_INT128_T
#endif

#define MD_XTC_CACHE_MAGIC   0x8281237612371
#define MD_XTC_CACHE_VERSION 3

#define MD_XTC_TRAJ_MAGIC 0x162365dac721995

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

typedef struct xtc_t {
    uint64_t magic;
    md_file_o* file;
    str_t filepath; // Store path to create multiple file streams
    md_array(int64_t) frame_offsets;
    md_trajectory_header_t header;
    md_allocator_i* allocator;
} xtc_t;

// XDR Specific stuff
static inline bool xdr_read_bytes (uint8_t* ptr, size_t count, md_file_o* xdr_file) {
    ASSERT(ptr);
    ASSERT(xdr_file);

    size_t bytes = count;
    if (md_file_read(xdr_file, ptr, bytes) == bytes) {
        return true;
    }

    return false;
}

static inline bool xdr_read_int32(int32_t* ptr, size_t count, md_file_o* xdr_file) {
    ASSERT(ptr);
    ASSERT(xdr_file);

    size_t bytes = count * sizeof(int32_t);
    if (md_file_read(xdr_file, ptr, bytes) == bytes) {
#if __LITTLE_ENDIAN__
        for (size_t i = 0; i < count; ++i) {
            ptr[i] = BSWAP32(ptr[i]);
        }
#endif
        return true;
    }

    return false;
}

static inline bool xdr_read_float (float* ptr, size_t count, md_file_o* xdr_file) {
    ASSERT(ptr);
    ASSERT(xdr_file);
    return xdr_read_int32((int32_t*)ptr, count, xdr_file);
}

static inline int sizeofint(int size) {
    unsigned int num = 1;
    int num_of_bits = 0;

    while (size >= (int)num && num_of_bits < 32) {
        num_of_bits++;
        num <<= 1;
    }
    return num_of_bits;
}

static inline int sizeofints(int num_of_ints, unsigned int sizes[]) {
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

static const uint32_t magicints[] = {
    0,        0,        0,       0,       0,       0,       0,       0,       0,       8,
    10,       12,       16,      20,      25,      32,      40,      50,      64,      80,
    101,      128,      161,     203,     256,     322,     406,     512,     645,     812,
    1024,     1290,     1625,    2048,    2580,    3250,    4096,    5060,    6501,    8192,
    10321,    13003,    16384,   20642,   26007,   32768,   41285,   52015,   65536,   82570,
    104031,   131072,   165140,  208063,  262144,  330280,  416127,  524287,  660561,  832255,
    1048576,  1321122,  1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042, 8388607,
    10568983, 13316085, 16777216};

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(magicints) / sizeof(*magicints))

#define BRANCHFREE_DIV

#ifdef BRANCHFREE_DIV
#define DIV_T struct libdivide_u64_branchfree_t
#define DIV_INIT(x) libdivide_u64_branchfree_gen(x)
#define DIV(x, y) libdivide_u64_branchfree_do(x, y)

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
#else
#define DIV_T struct libdivide_u64_t
#define DIV_INIT(x) libdivide_u64_gen(x)
#define DIV(x, y) libdivide_u64_do(x, y)

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

static const DIV_T denoms_64_2[] = {
    {0, 6}, {0X47AE147AE147AE15, 70}, {0XE38E38E38E38E38F, 7}, {0, 8}, {0X47AE147AE147AE15, 72}, {0XD1B71758E219652C, 9}, {0, 10}, {0X47AE147AE147AE15, 74}, {0XD1B71758E219652C, 11},
    {0, 12}, {0X47AE147AE147AE15, 76}, {0XCD953E4FF45B0BF9, 13}, {0, 14}, {0XA1CFA066147BB066, 14}, {0XCB90017DAE02CBA7, 15}, {0, 16}, {0XA1CFA066147BB066, 16}, {0XCB90017DAE02CBA7, 17},
    {0, 18}, {0XA14F47E0CDE8FB79, 18}, {0XCB90017DAE02CBA7, 19}, {0, 20}, {0XA14F47E0CDE8FB79, 20}, {0XCB4FE340A8C89403, 21}, {0, 22}, {0XA14F47E0CDE8FB79, 22}, {0XCB4FE340A8C89403, 23},
    {0, 24}, {0X4F7F449D26A949DD, 88}, {0XCB3FE06D28DF9EEF, 25}, {0, 26}, {0XA147476AF31DC745, 26}, {0XCB37DFB8F7D92618, 27}, {0, 28}, {0XA147476AF31DC745, 28}, {0XCB33DF8C3EE5EA77, 29},
    {0, 30}, {0X428A8ECA9A8C6397, 94}, {0XCB31DF8139CAE935, 31}, {0, 32}, {0X428A8ECA9A8C6397, 96}, {0XCB30DF7E8D046408, 33}, {0, 34}, {0X428A8ECA9A8C6397, 98}, {0XCB305F7DEC10CAE5, 35},
    {0, 36}, {0X428A8ECA9A8C6397, 100}, {0X96603EFB91E54C07, 101}, {0X40000C000201, 101}, {0X428A4ECA87BD4923, 102}, {0X965FFEFB8574EA53, 103}, {0, 40}, {0X428A4ECA87BD4923, 104},
    {0X965FFEFB8574EA53, 105}, {0, 42}, {0X428A3ECA8603773D, 106}, {0XCB2FF77DC25ACAB6, 43}, {0, 44}, {0XA1451B6542CC6CA8, 44}, {0XCB2FF77DC25ACAB6, 45}, {0X400000C0001, 109},
    {0XA145196542C00E8D, 46}, {0XCB2FF57DC255C633, 47}, {0, 48},
};
#endif


// Bitreader
typedef struct br_t {
    const uint64_t* stream;
    uint64_t data;
    uint64_t next;
    uint32_t cache_bits;
    uint32_t stream_size;
} br_t;

static inline void br_init(br_t* r, const uint64_t* stream, size_t num_qwords) {
    ASSERT(num_qwords >= 2);
    r->stream = stream + 2;
    r->data   = stream[0];
    r->next   = stream[1];
#if __LITTLE_ENDIAN__
    r->data   = BSWAP64(r->data);
    r->next   = BSWAP64(r->next);
#endif
    r->cache_bits = 128;
    r->stream_size = (uint32_t)num_qwords - 2;
}

static inline void br_load_next(br_t* r) {
    if (r->cache_bits > 64 || r->stream_size == 0) return;

    uint64_t data = *r->stream++;
    r->stream_size -= 1;

#if __LITTLE_ENDIAN__
    data = BSWAP64(data);
#endif

    // Always perform unified update
    uint64_t value = (r->cache_bits < 64) ? (data >> r->cache_bits) : 0;
    r->data |= value;
    r->next  = data << (64 - r->cache_bits);
    r->cache_bits += 64;
}

static inline uint64_t br_read(br_t* r, size_t num_bits) {
    ASSERT(num_bits <= 64);

    uint64_t shift = 64 - num_bits;
    uint64_t res   = r->data >> shift;
    r->cache_bits -= (uint32_t)num_bits;

    // Unified handling
    if (num_bits < 64) {
        r->data <<= num_bits;
        r->data |= r->next >> shift;
        r->next <<= num_bits;
    } else {
        r->data = r->next;
    }

    br_load_next(r);
    return res;
}

static inline uint64_t br_peek(br_t* r, size_t num_bits) {
    ASSERT(num_bits <= 64);
    uint64_t shift = 64 - num_bits;
    uint64_t res   = r->data >> shift;
    return res;
}

static inline void br_skip(br_t* r, size_t num_bits) {
    ASSERT(num_bits < 64);
    r->cache_bits -= (uint32_t)num_bits;
    r->data <<= num_bits;

    // Append extracted bits from next
    r->data |= r->next >> (64 - num_bits);
    r->next <<= num_bits;

    br_load_next(r);
}

typedef struct {
    uint8_t  num_of_bits;
    uint8_t  part_mask;
    uint8_t  big_shift;
    uint8_t  sml_shift;
} bit_data_t;

typedef struct {
    uint32_t size_y;
    uint32_t size_z;
    bit_data_t bit;
    DIV_T    div_zy;
    DIV_T    div_z;
} unpack_data_t;

typedef md_128i v4i_t;
#define v4i_set(x, y, z, w) md_mm_set_epi32(w, z, y, x)
#define v4i_set1(x)         md_mm_set1_epi32(x)
#define v4i_add(a, b)       md_mm_add_epi32(a, b)
#define v4i_sub(a, b)       md_mm_sub_epi32(a, b)
#define v4i_load(src)       md_mm_loadu_epi32(src)

static inline void write_coord(float* dst, v4i_t coord, md_128 invp) {
    md_128  data = md_mm_mul_ps(md_mm_cvtepi32_ps(coord), invp);
    MEMCPY(dst, &data, sizeof(float) * 3);
}

static inline void init_unpack_bit_data(bit_data_t* data, uint32_t num_of_bits) {
    uint32_t partbits = num_of_bits & 7;
    data->num_of_bits = (uint8_t)num_of_bits;
    data->big_shift = 64 - ((num_of_bits + 7) & ~7);  // Align to the next multiple of 8
    data->sml_shift = (8 - partbits) & 7;
    data->part_mask = partbits ? (1ul << partbits) - 1 : 0xFF;
}

static inline md_128i md_mm_bswap_128(md_128i v) {
    const md_128i shuf_mask = md_mm_set_epi8(
        0x0, 0x1, 0x2, 0x3,
        0x4, 0x5, 0x6, 0x7,
        0x8, 0x9, 0xA, 0xB,
        0xC, 0xD, 0xE, 0xF
    );
    md_128i swapped = simde_mm_shuffle_epi8(v, shuf_mask);
    return swapped;
}

static inline md_128i md_mm_bswap_64(md_128i v) {
    const md_128i shuf_mask = md_mm_set_epi8(
        0x7, 0x6, 0x5, 0x4,
        0x3, 0x2, 0x1, 0x0,
        0xF, 0xE, 0xD, 0xC,
        0xB, 0xA, 0x9, 0x8
    );
    md_128i swapped = simde_mm_shuffle_epi8(v, shuf_mask);
    return swapped;
}

// THIS ASSUMES BIG-ENDIAN INPUT DATA
// AND BIT_LENGTH IS > 64 AND <= 121
static inline void extract_bits_be_raw_121(uint64_t out[2],
    const uint8_t* base,
    size_t bit_offset,
    size_t bit_length)   // 1..121
{
    ASSERT(64 < bit_length && bit_length <= 121);

    size_t byte_offset = bit_offset >> 3;
    size_t bit_in_byte = bit_offset & 7;
    int shift_right = (int)(128 - bit_length);

#if defined(MD_XTC_VECTOR_EXTRACT)
    md_128i v = md_mm_loadu_si128((const md_128i*)(base + byte_offset));
    v = md_mm_bswap_64(v);

    md_128i shift = md_mm_set1_epi64((int64_t)bit_in_byte);
    md_128i shifted = simde_mm_sll_epi64(v, shift);

    md_128i carry = simde_mm_srl_epi64(v, md_mm_set1_epi64(64 - (int64_t)bit_in_byte));
    carry = simde_mm_slli_si128(carry, 8);

    md_128i combined = simde_mm_or_si128(shifted, carry);

    uint64_t out[2];
    md_mm_storeu_si128((md_128i*)out, combined);

    out[1] = out[1] >> shift_right;
#else
    uint64_t raw[2];
    MEMCPY(raw, base + byte_offset, sizeof(raw));
    uint64_t swapped0 = BSWAP64(raw[0]);
    uint64_t swapped1 = BSWAP64(raw[1]);
    uint64_t high = (swapped0 << bit_in_byte);
    uint64_t low  = (swapped1 << bit_in_byte) | (swapped0 >> (64 - bit_in_byte));
    out[0] = high;
    out[1] = low >> shift_right;
#endif
}

// Loads 16 bytes from base+(bit_offset>>3), byteswaps, then left-shifts the 128-bit window by (bit_offset&7).
// out[0] = high 64 bits, out[1] = low 64 bits, both MSB-first aligned to bit_offset.
static inline void load_shifted_be128(uint64_t out[2], const uint8_t* base, size_t bit_offset) {
    const size_t byte_offset = bit_offset >> 3;
    const size_t bit_in_byte = bit_offset & 7;

    uint64_t raw[2];
    MEMCPY(raw, base + byte_offset, sizeof(raw));
    uint64_t hi = BSWAP64(raw[0]);
    uint64_t lo = BSWAP64(raw[1]);

    const uint64_t k = bit_in_byte;
    const uint64_t nz_mask = -(uint64_t)(k != 0);
    const uint64_t s = (64u - (uint64_t)k) & 63u; // never 64

    uint64_t high = (hi << k) | ((lo >> s) & nz_mask);
    uint64_t low  = (lo << k);

    out[0] = high;
    out[1] = low;
}

// Extract up to 32 bits starting at `start` (0..127) where 0 is MSB of out[0].
static inline uint32_t get_bits_be128_u32(const uint64_t win[2], uint32_t start, uint32_t len) {
    ASSERT(len > 0 && len <= 32);
    if (start < 64) {
        const uint64_t x = win[0] << start;
        if (start + len <= 64) {
            return (uint32_t)(x >> (64 - len));
        } else {
            // crosses into win[1]
            const uint32_t a = 64 - start;
            const uint32_t b = len - a;
            const uint32_t hi_part = (uint32_t)(x >> (64 - a));
            const uint32_t lo_part = (uint32_t)(win[1] >> (64 - b));
            return (hi_part << b) | lo_part;
        }
    } else {
        const uint32_t s = start - 64;
        return (uint32_t)((win[1] << s) >> (64 - len));
    }
}

static inline uint64_t extract_bits_be_raw_64(
    const uint8_t* base,
    size_t bit_offset,
    size_t bit_length)   // 1..64
{
    size_t byte_offset = bit_offset >> 3;
    size_t bit_in_byte = bit_offset & 7;
    
    // Always read 16 bytes (two uint64_t)
    uint64_t raw[2];
    MEMCPY(raw, base + byte_offset, sizeof(raw));
    
    uint64_t hi = BSWAP64(raw[0]);
    uint64_t lo = BSWAP64(raw[1]);
    
    // Combine (lo part masked out if not needed)
    uint64_t nz_mask = -(uint64_t)(bit_in_byte != 0);
    uint64_t combined = (hi << bit_in_byte) | ((lo >> (64 - bit_in_byte)) & nz_mask);
    
    return combined >> (64 - bit_length);
}

static inline uint64_t extract_bits_be_raw_57(
    const uint8_t* base,
    size_t bit_offset,
    size_t bit_length)   // 1..64
{
    // Calculate byte offset
    size_t byte_offset = bit_offset >> 3;
    size_t bit_in_byte = bit_offset & 7;

    // Read 8 bytes unaligned (covers up to 64 bits + 7 bit offset)
    uint64_t raw;
    MEMCPY(&raw, base + byte_offset, sizeof(raw));

    // Byteswap to host endian (assuming host is little-endian)
    uint64_t swapped = BSWAP64(raw);

    // Extract the bits (they're now at the top after bswap)
	// Mask upper by shifting left and then right
    return (swapped << bit_in_byte) >> (64 - bit_length);
}

static inline uint32_t extract_bits_be_raw_32(
    const uint8_t* base,
    size_t bit_offset,
    size_t bit_length)   // 1..32
{
    size_t byte_offset = bit_offset >> 3;
    size_t bit_in_byte = bit_offset & 7;
    
    // Read 8 bytes (or could read just 5 bytes for 32-bit case)
    uint64_t raw;
    MEMCPY(&raw, base + byte_offset, sizeof(raw));
    
    uint64_t swapped = BSWAP64(raw);
    return (uint32_t)((swapped << bit_in_byte) >> (64 - bit_length));
}

// Method for 25 bits extraction
static inline uint32_t extract_bits_be_raw_25(
    const uint8_t* base,
    size_t bit_offset,
    size_t bit_length)   // 1..32
{
    size_t byte_offset = bit_offset >> 3;
    size_t bit_in_byte = bit_offset & 7;

    // Read 4 bytes unaligned
    uint32_t raw;
    MEMCPY(&raw, base + byte_offset, sizeof(raw));

    uint32_t swapped = BSWAP32(raw);
    return (uint32_t)((swapped << bit_in_byte) >> (32 - bit_length));
}

static inline v4i_t unpack_coord64(uint64_t w, const unpack_data_t* unpack) {
    uint64_t t = ((w << unpack->bit.sml_shift) & (~0xFFull)) | (w & unpack->bit.part_mask);
    uint64_t v = BSWAP64(t) >> unpack->bit.big_shift;

    uint32_t x  = DIV(v, &unpack->div_zy);   // v / (Y*Z)
    uint64_t yz = DIV(v, &unpack->div_z);    // v / Z

    uint32_t y = yz - (uint64_t)x  * unpack->size_y;
    uint32_t z = v  - (uint64_t)yz * unpack->size_z;

	return v4i_set(x, y, z, 0);
}


static inline v4i_t unpack_coord128(uint64_t w[2], const unpack_data_t* unpack) {
    int fullbytes = unpack->bit.num_of_bits >> 3;
    int partbits  = unpack->bit.num_of_bits  & 7;
    const uint64_t zy = (uint64_t)unpack->size_z * unpack->size_y;

#ifdef HAS_INT128_T
    __uint128_t v;
    MEMCPY(&v, w, sizeof(v));
    uint32_t x = (uint32_t)(v / zy);
    uint64_t q = (uint64_t)(v - x*zy);
    uint32_t y = (uint32_t)(q / unpack->size_z);
    uint32_t z = (uint32_t)(q % unpack->size_z);
    return v4i_set(x, y, z, 0);
#elif defined(__x86_64__) && defined(_MSC_VER) && (_MSC_VER >= 1920)
    uint64_t q = 0;
    uint32_t x = (uint32_t)_udiv128(w[1], w[0], zy, &q);
    uint32_t y = (uint32_t)(q / unpack->size_z);
    uint32_t z = (uint32_t)(q % unpack->size_z);
    return v4i_set(x, y, z, 0);
#else
    // Default fallback
    // limbs[0] is least-significant 32 bits
    const int num_of_bytes = fullbytes + (partbits ? 1 : 0);
    const uint32_t sizes[3] = {0, unpack->size_y, unpack->size_z};
    uint32_t nums[4]  = {0};
    uint32_t limbs[4] = {0};
    MEMCPY(limbs, w, (size_t)num_of_bytes);

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

static inline v4i_t read_ints3(const uint8_t* base, size_t bit_offset, const bit_data_t bit_data[3]) {
    uint32_t val[4];
    int offset = 0;
    for (int i = 0; i < 3; ++i) {
        int num_bits  = bit_data[i].num_of_bits;
        int big_shift = bit_data[i].big_shift;
        int sml_shift = bit_data[i].sml_shift;
        int part_mask = bit_data[i].part_mask;
        uint32_t w = extract_bits_be_raw_32(base, bit_offset + offset, num_bits);
        uint32_t r = ((w << sml_shift) & (~0xFFul)) | (w & part_mask);
        uint32_t v = BSWAP32(r) >> big_shift;
        val[i] = v;
        offset += num_bits;
    }
    return v4i_load(val);
}

bool md_xtc_read_frame_header(md_file_o* xdr, int* natoms, int* step, float* time, float box[3][3]) {
    int magic;

    if (!xdr_read_int32(&magic, 1, xdr)) {
        MD_LOG_ERROR("XTC: Failed to read magic number in header");
        return false;
    }
    if (magic != XTC_MAGIC) {
        MD_LOG_ERROR("XTC: Magic number did not match");
        return false;
    }
    if (!xdr_read_int32(natoms, 1, xdr)) {
        MD_LOG_ERROR("XTC: Failed to read number of atoms");
        return false;
    }
    if (!xdr_read_int32(step, 1, xdr)) {
        MD_LOG_ERROR("XTC: Failed to read step");
        return false;
    }
    if (!xdr_read_float(time, 1, xdr)) {
        MD_LOG_ERROR("XTC: Failed to read timestamp");
        return false;
    }
    if (!xdr_read_float((float*)box, 9, xdr)) {
        MD_LOG_ERROR("XTC: Failed to read box dimensions");
        return false;
    }

    return true;
}

size_t md_xtc_read_frame_offsets_and_times(md_file_o* xdr, md_array(int64_t)* frame_offsets, md_array(double)* frame_times, md_allocator_i* alloc) {
    int step, natoms;
    float time;
    float box[3][3];

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

    if (!md_xtc_read_frame_header(xdr, &natoms, &step, &time, box)) {
        MD_LOG_ERROR("XTC: File does not appear to be a valid xtc trajectory");
        return 0;
    }

    size_t num_frames = 0;

    /* Dont bother with compression for nine atoms or less */
    if (natoms <= 9) {
        const size_t framebytes = XTC_SMALL_HEADER_SIZE + XTC_SMALL_COORDS_SIZE * natoms;
        const size_t est_frames = (filesize / framebytes); /* Should we complain if framesize doesn't divide filesize? */

        md_array_ensure(*frame_offsets, est_frames, alloc);
        md_array_ensure(*frame_times,   est_frames, alloc);

        // Push first frame
        md_array_push(*frame_offsets, 0, alloc);
        md_array_push(*frame_times, time, alloc);
        num_frames += 1;

        for (size_t i = 1; i < est_frames; i++) {
            const size_t offset = i * framebytes;

            if (!md_file_seek(xdr, offset, MD_FILE_BEG) || !md_xtc_read_frame_header(xdr, &natoms, &step, &time, box)) {
                MD_LOG_DEBUG("XTC: encountered corrupted frame header");
                break;
            }

            // Push frame i
            md_array_push(*frame_offsets, offset, alloc);
            md_array_push(*frame_times, time, alloc);
            num_frames += 1;
        }
        md_array_push(*frame_offsets, num_frames * framebytes, alloc);
    } else {
        int framebytes, est_nframes;

        /* Move pos back to end of first header */
        if (!md_file_seek(xdr, XTC_HEADER_SIZE, MD_FILE_BEG)) {
            return 0;
        }

        if (!xdr_read_int32(&framebytes, 1, xdr)) {
            MD_LOG_ERROR("XTC: Failed to read framebytes");
            return 0;
        }
        framebytes = (framebytes + 3) & ~0x03; /* Rounding to the next 32-bit boundary */
        est_nframes = (int)(filesize / ((int64_t)(framebytes + XTC_HEADER_SIZE)) + 1);
        /* First `framebytes` might be larger than average, so we would underestimate `est_nframes`*/
        est_nframes += est_nframes / 5;

        /* Skip `framebytes` */
        if (!md_file_seek(xdr, (int64_t)(framebytes), MD_FILE_CUR)) {
            MD_LOG_DEBUG("XTC: encountered corrupted frame");
        }

        md_array_ensure(*frame_offsets, (size_t)est_nframes, alloc);
        md_array_ensure(*frame_times,   (size_t)est_nframes, alloc);

        md_array_push(*frame_offsets, 0, alloc);
        md_array_push(*frame_times, time, alloc);
        num_frames += 1;

        while (true) {
            const int64_t offset = md_file_tell(xdr);
            if (offset == (int64_t)filesize) {
                break;
            }
            if (!md_xtc_read_frame_header(xdr, &natoms, &step, &time, box)) {
                MD_LOG_DEBUG("XTC: encountered corrupted frame header");
                goto done;
            }

            // Skip natoms + additional frame data
            md_file_seek(xdr, 4 + 32, MD_FILE_CUR);

            /* Read how much to skip */
            if (!xdr_read_int32(&framebytes, 1, xdr)) {
                MD_LOG_DEBUG("XTC: encountered corrupted frame");
                goto done;
            }
            framebytes = (framebytes + 3) & ~0x03; /* Rounding to the next 32-bit boundary */

            /* Skip `framebytes` to next header */
            if (!md_file_seek(xdr, framebytes, MD_FILE_CUR)) {
                MD_LOG_DEBUG("XTC: encountered corrupted frame");
                goto done;
            }

            /* Store position in `offsets`, adjust for header */
            md_array_push(*frame_offsets, offset, alloc);
            md_array_push(*frame_times, time, alloc);
            num_frames += 1;
        }
        // Add last offset
        md_array_push(*frame_offsets, filesize, alloc);
    }
done:
    return num_frames;
}

//#define PRINT_DEBUG

#define BUFFER_SIZE KILOBYTES(64)
#define ALIGNMENT 64
#define GUARD 16

#if MD_PLATFORM_UNIX
#include <stdio.h>
#include <fcntl.h>
#endif

#if MD_PLATFORM_WINDOWS
#define ALIGNED_ALLOC(size, alignment) _aligned_malloc(size, alignment)
#define ALIGNED_FREE(ptr) _aligned_free(ptr)
#else
#define ALIGNED_ALLOC(size, alignment) aligned_alloc(alignment, size)
#define ALIGNED_FREE(ptr) free(ptr)
#endif

typedef struct bytestream_t {
    md_file_o* file;
    uint8_t* buf;
    size_t pos;   // next byte to consume
    size_t end;   // one past last valid byte
} bytestream_t;

static inline bool bytestream_init(bytestream_t* bs, md_file_o* file) {
    bs->file = file;
    bs->pos = 0;
    bs->end = 0;
    bs->buf = ALIGNED_ALLOC(BUFFER_SIZE, ALIGNMENT);
    if (bs->buf == NULL) {
        return false;
    }
    return true;
}

static inline void bytestream_free(bytestream_t* bs) {
    if (bs->buf) {
        ALIGNED_FREE(bs->buf);
        bs->buf = NULL;
    }
}

static inline int bytestream_refill(bytestream_t *s) {
    // Move remaining data to front
    size_t remain = s->end - s->pos;
    MEMMOVE(s->buf, s->buf + s->pos, remain);

    s->pos = 0;
    s->end = remain;
	size_t n = md_file_read(s->file, s->buf + remain, BUFFER_SIZE - remain);
    if (n == 0) return n;

    s->end += n;
    return 1;
}

static inline uint8_t *bytestream_ensure(bytestream_t *s, size_t need) {
    while (s->end - s->pos < need) {
        if (!bytestream_refill(s)) break;
    }
    return s->buf + s->pos;
}

static inline void bytestream_advance(bytestream_t *s, size_t n) {
    s->pos += n;
}

static inline void read_ints(int* dst, size_t count, uint8_t* src) {
	MEMCPY(dst, src, count * sizeof(int32_t));
#if __LITTLE_ENDIAN__
    for (size_t i = 0; i < count; i++) {
        dst[i] = BSWAP32(dst[i]);
    }
#endif
}

#define USE_BYTESTREAM 0

size_t md_xtc_read_frame_coords(md_file_o* xdr_file, float* out_coords_ptr, size_t out_coord_cap) {
    if (xdr_file == NULL || out_coords_ptr == NULL) {
        return 0;
    }
    
#if USE_BYTESTREAM
	bytestream_t stream = { 0 };
    if (!bytestream_init(&stream, xdr_file)) {
        MD_LOG_ERROR("XTC: Failed to allocate bytestream buffer");
        return 0;
	}
    uint8_t* base_ptr = bytestream_ensure(&stream, 64);
#endif

    int atom_idx = 0;

    int natoms;
#if USE_BYTESTREAM
    read_ints(&natoms, 1, base_ptr);
    base_ptr += 4;
#else
    if (!xdr_read_int32(&natoms, 1, xdr_file)) {
        return 0;
    }
#endif

    if (out_coord_cap < (size_t)natoms) {
        MD_LOG_ERROR("The supplied coordinate buffer is not sufficient to fit the number of atoms within the frame");
        return 0;
    }

    /* Dont bother with compression for three atoms or less */
    if (natoms <= 9) {
#if USE_BYTESTREAM
		read_ints((int*)out_coords_ptr, natoms * 3, base_ptr);
#else
        if (!xdr_read_float(out_coords_ptr, natoms * 3, xdr_file)) {
            return 0;
        }
#endif
        return (size_t)natoms;
    }

    /* Compression-time if we got here. Read precision first */
    float precision;
#if USE_BYTESTREAM
    read_ints((int*)&precision, 1, base_ptr);
	base_ptr += 4;
#else
    if (!xdr_read_float(&precision, 1, xdr_file)) {
        return 0;
    }
#endif

    int32_t minint[3], maxint[3];
#if USE_BYTESTREAM
    read_ints(minint, 3, base_ptr);
    base_ptr += 12;
    read_ints(maxint, 3, base_ptr);
    base_ptr += 12;
#else
    if (!xdr_read_int32(minint, 3, xdr_file) || !xdr_read_int32(maxint, 3, xdr_file)) {
        return 0;
    }
#endif

    uint32_t sizeint[3] = {
        (uint32_t)(maxint[0] - minint[0] + 1),
        (uint32_t)(maxint[1] - minint[1] + 1),
        (uint32_t)(maxint[2] - minint[2] + 1),
    };

    uint32_t bitsize = 0;
    uint32_t big_bits = 0;
    unpack_data_t big_unpack = {0};
    bit_data_t big_bit_data[3] = {0};

    /* check if one of the sizes is to big to be multiplied */
    if ((sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffff) {
        for (int i = 0; i < 3; ++i) {
            int bits = sizeofint(sizeint[i]);
            init_unpack_bit_data(&big_bit_data[i], bits);
            big_bits += bits;
        }
    } else {
        bitsize = sizeofints(3, sizeint);
        big_bits = bitsize;
        big_unpack.size_y = sizeint[1];
        big_unpack.size_z = sizeint[2];
        big_unpack.div_zy = DIV_INIT((uint64_t)sizeint[1] * sizeint[2]);
        big_unpack.div_z  = DIV_INIT(sizeint[2]);
        init_unpack_bit_data(&big_unpack.bit, bitsize);
    }

    int smallidx;
#if USE_BYTESTREAM
    read_ints(&smallidx, 1, base_ptr);
	base_ptr += 4;
#else
    if (!xdr_read_int32(&smallidx, 1, xdr_file)) {
        return 0;
    }
#endif

    int smaller = (smallidx > FIRSTIDX) ? (int)magicints[smallidx - 1] / 2 : 0;
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
    uint32_t num_bytes = 0;
#if USE_BYTESTREAM
    read_ints((int*)&num_bytes, 1, base_ptr);
    base_ptr += 4;
#else
    if (!xdr_read_int32((int32_t*)&num_bytes, 1, xdr_file)) {
        return 0;
    }
    size_t mem_size = ALIGN_TO(num_bytes, 32);
    char* mem = malloc(mem_size);
    uint8_t* base = mem;

    if (!xdr_read_bytes((uint8_t*)mem, num_bytes, xdr_file)) {
        goto done;
    }
#endif

    float* lfp = out_coords_ptr;
	const md_128 invp = md_mm_set1_ps(1.0f / precision);
    const v4i_t vminint = v4i_set(minint[0], minint[1], minint[2], 0);
    v4i_t thiscoord;
    int run = 0;
    int run_count = 0;

#ifdef PRINT_DEBUG
    uint32_t smlbits_hist[LASTIDX] = {0};
    uint32_t rle_hist[9] = {0};
    MD_LOG_DEBUG("bigsize: %i", bitsize);
#endif

#if USE_BYTESTREAM
	// Advance to compressed base data
	bytestream_advance(&stream, 10 * sizeof(int32_t));
#endif
    size_t bit_offset = 0;

    while (atom_idx < natoms) {
#if USE_BYTESTREAM
        uint8_t* base = bytestream_ensure(&stream, 256);
#endif
        uint32_t run_data = extract_bits_be_raw_25(base, bit_offset + big_bits, 6);
        uint64_t raw[2];
        load_shifted_be128(raw, base, bit_offset);

        if (bitsize == 0) {
            thiscoord = read_ints3(base, bit_offset, big_bit_data);
        } else {
            if (big_bits <= 64) {
                thiscoord = unpack_coord64(raw[0], &big_unpack);
            } else {
                thiscoord = unpack_coord128(raw, &big_unpack);
            }
        }

        thiscoord = v4i_add(thiscoord, vminint);

        uint32_t run_flag = run_data & 32;
        uint32_t run_skip = run_flag ? 6 : 1;
        bit_offset += big_bits + run_skip;

        int is_smaller = 0;
        if (run_flag) {
            run = run_data & 31;
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

#ifdef PRINT_DEBUG
        rle_hist[run_count]++;
#endif

        // If run, write first atom after second atom
        int offset = run ? 3 : 0;
        write_coord(lfp + offset, thiscoord, invp); lfp += 3;

        if (run) {
            const v4i_t vsmall = v4i_set1(smallnum);
            const int sml_bits = sml_unpack.bit.num_of_bits;

            uint64_t raw;
            if (sml_bits <= 57) {
                raw = extract_bits_be_raw_57(base, bit_offset, sml_bits);
            } else {
                raw = extract_bits_be_raw_64(base, bit_offset, sml_bits);
			}
            
            v4i_t coord = unpack_coord64(raw, &sml_unpack);
            thiscoord = v4i_add(coord, v4i_sub(thiscoord, vsmall));
            
            // Write second atom before first atom
            write_coord(lfp - 3, thiscoord, invp); lfp += 3;
            bit_offset += sml_bits;

            for (int i = 1; i < run_count; i += 1) {
                if (sml_bits <= 57) {
                    raw = extract_bits_be_raw_57(base, bit_offset, sml_bits);
                } else {
                    raw = extract_bits_be_raw_64(base, bit_offset, sml_bits);
                }
                coord = unpack_coord64(raw, &sml_unpack);
                thiscoord = v4i_add(coord, v4i_sub(thiscoord, vsmall));
                write_coord(lfp, thiscoord, invp); lfp += 3;

                bit_offset += sml_bits;
            }
        }

        if (is_smaller) {
            smallidx += is_smaller;
            smallnum = (int)magicints[smallidx] / 2;
            smaller  = (smallidx > FIRSTIDX) ? (int)magicints[smallidx - 1] / 2 : 0;
        }
        if (smallidx != sml_unpack.bit.num_of_bits) {
#ifdef PRINT_DEBUG
            smlbits_hist[smallidx]++;
#endif
            uint32_t sml_size       = magicints[smallidx];
            sml_unpack.size_y       = sml_size;
            sml_unpack.size_z       = sml_size;
            sml_unpack.div_zy       = denoms_64_2[smallidx - FIRSTIDX];
            sml_unpack.div_z        = denoms_64_1[smallidx - FIRSTIDX];
            init_unpack_bit_data(&sml_unpack.bit, smallidx);
        }
        if (smallidx < FIRSTIDX) {
            MD_LOG_ERROR("XTC: Invalid size found in 'xdrfile_decompress_coord_float'.");
            goto done;
        }

#if USE_BYTESTREAM
		size_t advanced_bytes = bit_offset / 8;
		bytestream_advance(&stream, advanced_bytes);
		bit_offset &= 7;
#endif
    }
    
done:
#if USE_BYTESTREAM
	bytestream_free(&stream);
#else
    free(mem);
#endif

#ifdef PRINT_DEBUG
    int max_smallbits = 0;
    uint32_t total_smallbits = 0;
    for (int i = 0; i < LASTIDX; i++) {
        if (smlbits_hist[i] > smlbits_hist[max_smallbits]) {
            max_smallbits = i;
        }
        total_smallbits += smlbits_hist[i];
    }

    uint32_t total_rle = 0;
    for (int i = 0; i < 9; ++i) {
        total_rle += rle_hist[i];
    }

    // Printf a histogram with horizontal bars
    printf("--- Smallbits ---\n");
    for (int i = 0; i < LASTIDX; i++) {
        if (smlbits_hist[i] > 0) {
            int num_hashes = (int)(smlbits_hist[i] * 50 / total_smallbits);
            printf("%2d: ", i);
            for (int j = 0; j < num_hashes; j++) {
                printf("#");
            }
            printf(" (%d)\n", smlbits_hist[i]);
        }
    }
    printf("----------------\n");
    printf("--- run length ---\n");
    for (int i = 0; i < 9; i++) {
        if (rle_hist[i] > 0) {
            int num_hashes = (int)(rle_hist[i] * 50 / total_rle);
            printf("%2d: ", i);
            for (int j = 0; j < num_hashes; j++) {
                printf("#");
            }
            printf(" (%d)\n", rle_hist[i]);
        }
    }
    printf("----------------\n\n");
#endif

    return (size_t)atom_idx;
}

bool xtc_get_header(struct md_trajectory_o* inst, md_trajectory_header_t* header) {
    xtc_t* xtc = (xtc_t*)inst;
    ASSERT(xtc);
    ASSERT(xtc->magic == MD_XTC_TRAJ_MAGIC);
    ASSERT(header);

    *header = xtc->header;
    return true;
}

static bool xtc_decode_frame_data(md_file_o* file, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(file);

    bool result = true;

    if ((x || y || z) && !(x && y && z)) {
        MD_LOG_ERROR("XTC: User supplied coordinates (x,y,z) cannot be partially supplied");
        return false;
    }

    // Get header
    int natoms = 0, step = 0;
    float time = 0;
    float box[3][3];
    result = md_xtc_read_frame_header(file, &natoms, &step, &time, box);
    if (result) {
        if (header) {
            // nm -> Ångström
            for (int i = 0; i < 3; ++i) {
                box[i][0] *= 10.0f;
                box[i][1] *= 10.0f;
                box[i][2] *= 10.0f;
            }
            header->num_atoms = natoms;
            header->index = step;
            header->timestamp = time;
            header->unitcell = md_unitcell_from_matrix_float(box);
        }

        if (x && y && z) {
            size_t byte_size = natoms * sizeof(float) * 3;
            float* xyz = md_alloc(md_get_heap_allocator(), byte_size);
            size_t written_count = md_xtc_read_frame_coords(file, xyz, natoms);
            result = (written_count == natoms);
            if (result) {
                // nm -> Ångström
                for (int i = 0; i < natoms; ++i) {
                    x[i] = xyz[i*3 + 0] * 10.0f;
                    y[i] = xyz[i*3 + 1] * 10.0f;
                    z[i] = xyz[i*3 + 2] * 10.0f;
                }
            }
            md_free(md_get_heap_allocator(), xyz, byte_size);
        }
    }

    return result;
}

bool xtc_load_frame(struct md_trajectory_o* inst, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(inst);

    xtc_t* xtc = (xtc_t*)inst;
    if (xtc->magic != MD_XTC_TRAJ_MAGIC) {
        MD_LOG_ERROR("XTC: Error when decoding frame coord, xtc magic did not match");
        return false;
    }

    bool result = false;
    if (xtc->file) {
        md_file_o* file = md_file_open(xtc->filepath, MD_FILE_READ | MD_FILE_BINARY);
        if (!file) {
            MD_LOG_ERROR("XTC: Failed to open file");
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

        const int64_t pos = xtc->frame_offsets[frame_idx];
        if (!md_file_seek(file, pos, MD_FILE_BEG)) {
            MD_LOG_ERROR("XTC: Failed to seek to frame offset");
            md_file_close(file);
            return false;
        }

        result = xtc_decode_frame_data(file, header, x, y, z);
        md_file_close(file);
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
    md_file_o* file = md_file_open(cache_file, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
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
        md_file_close(file);
    }
    return result;
}

static bool write_cache(const xtc_cache_t* cache, str_t cache_file) {
    bool result = false;

    md_file_o* file = md_file_open(cache_file, MD_FILE_WRITE | MD_FILE_BINARY);
    if (!file) {
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
    md_file_close(file);
    return result;
}

md_trajectory_i* md_xtc_trajectory_create(str_t filename, md_allocator_i* ext_alloc, uint32_t flags) {
    ASSERT(ext_alloc);
    md_allocator_i* alloc = md_arena_allocator_create(ext_alloc, MEGABYTES(1));

    str_t path = str_copy(filename, alloc);
    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);

    if (file) {
        const size_t filesize = md_file_size(file);

        /* Go to file beg */
        if (!md_file_seek(file, 0, MD_FILE_BEG)) {
            MD_LOG_ERROR("XTC: Failed to seek to beginning of file");
            return 0;
        }

        int num_atoms, step;
        float time;
        float box[3][3];
        if (!md_xtc_read_frame_header(file, &num_atoms, &step, &time, box)) {
            goto fail;
        }

        if (num_atoms == 0) {
            MD_LOG_ERROR("XTC: Number of atoms in trajectory was zero");
            goto fail;
        }

        md_strb_t sb = md_strb_create(md_get_temp_allocator());
        md_strb_push_str(&sb, path);
        md_strb_push_str(&sb, STR_LIT(".cache"));
        str_t cache_path = md_strb_to_str(sb);

        xtc_cache_t cache = {0};
        if (!try_read_cache(&cache, cache_path, filesize, alloc)) {
            cache.header.magic     = MD_XTC_CACHE_MAGIC;
            cache.header.version   = MD_XTC_CACHE_VERSION;
            cache.header.num_bytes = filesize;
            cache.header.num_atoms = num_atoms;
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

        size_t max_frame_size = 0;
        for (size_t i = 0; i < cache.header.num_frames; ++i) {
            const size_t frame_size = (size_t)MAX(0, cache.frame_offsets[i + 1] - cache.frame_offsets[i]);
            max_frame_size = MAX(max_frame_size, frame_size);
        }

        void* mem = md_alloc(alloc, sizeof(md_trajectory_i) + sizeof(xtc_t));
        ASSERT(mem);
        MEMSET(mem, 0, sizeof(md_trajectory_i) + sizeof(xtc_t));

        md_trajectory_i* traj = mem;
        xtc_t* xtc = (xtc_t*)(traj + 1);

        xtc->magic = MD_XTC_TRAJ_MAGIC;
        xtc->allocator = alloc;
        xtc->file = file;
        xtc->filepath = path;
        xtc->frame_offsets = cache.frame_offsets;

        xtc->header = (md_trajectory_header_t) {
            .num_frames = cache.header.num_frames,
            .num_atoms = num_atoms,
            .max_frame_data_size = max_frame_size,
            .time_unit = md_unit_pikosecond(),
            .frame_times = cache.frame_times,
        };

        traj->inst = (struct md_trajectory_o*)xtc;
        traj->get_header = xtc_get_header;
        traj->load_frame = xtc_load_frame;

        return traj;
    }
fail:
    if (file) md_file_close(file);
    md_arena_allocator_destroy(alloc);
    return NULL;
}

void md_xtc_trajectory_free(md_trajectory_i* traj) {
    ASSERT(traj);
    ASSERT(traj->inst);
    xtc_t* xtc = (xtc_t*)traj->inst;
    if (xtc->magic != MD_XTC_TRAJ_MAGIC) {
        MD_LOG_ERROR("XTC: Cannot free trajectory, is not a valid XTC trajectory.");
        ASSERT(false);
        return;
    }

    if (xtc->file) md_file_close(xtc->file);
    md_arena_allocator_destroy(xtc->allocator);
}

static md_trajectory_loader_i xtc_loader = {
    md_xtc_trajectory_create,
    md_xtc_trajectory_free,
};

md_trajectory_loader_i* md_xtc_trajectory_loader(void) {
    return &xtc_loader;
}

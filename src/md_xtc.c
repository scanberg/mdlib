#include <md_xtc.h>

#include <md_util.h>
#include <md_trajectory.h>

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


// Bitreader
typedef struct br2_t {
    const uint8_t* base;
    size_t bit_offset;
} br2_t;

static inline uint64_t br2_peek_be(br2_t* r, size_t num_bits) {
    ASSERT(num_bits <= 64);
    size_t byte_offset = r->bit_offset >> 3;
    size_t shift = r->bit_offset & 7;
    uint64_t raw;
    MEMCPY(&raw, r->base + byte_offset, 8);
    uint64_t x = BSWAP64(raw);
    x <<= shift;
    x >>= 64 - num_bits;
    return x;
}

static inline uint64_t br2_read_be(br2_t* r, size_t num_bits) {
    uint64_t val = br2_peek_be(r, num_bits);
    r->bit_offset += num_bits;
    return val;
}

static inline void br2_skip(br2_t* r, size_t num_bits) {
    r->bit_offset += num_bits;
}

static inline uint64_t extract_u64_be(const uint8_t* base, size_t bit_offset, size_t num_bits) {
    ASSERT(num_bits <= 64);
    size_t byte_offset = bit_offset >> 3;
    size_t shift = bit_offset & 7;
    uint64_t x;
    MEMCPY(&x, base + byte_offset, 8);
    x = BSWAP64(x);
    x <<= shift;
    x >>= 64 - num_bits;
    return x;
}

typedef struct {
    uint64_t size_zy;
    uint32_t size_z;
    uint8_t  num_of_bits;
    uint8_t  part_mask;
    uint8_t  big_shift;
    uint8_t  sml_shift;
    uint64_t hi_shift;
    uint64_t lo_shift;
    uint64_t hi_mask;
    uint64_t lo_mask;
    DIV_T    div_zy;
    DIV_T    div_z;
} unpack_data_t;

typedef md_128i v4i_t;
#define v4i_set(x, y, z, w) md_mm_set_epi32(w, z, y, x)
#define v4i_set1(x)         md_mm_set1_epi32(x)
#define v4i_add(a, b)       md_mm_add_epi32(a, b)
#define v4i_sub(a, b)       md_mm_sub_epi32(a, b)
#define v4i_load(addr)      md_mm_loadu_epi32(addr)

static inline void write_coord(float* dst, v4i_t coord, float inv_precision) {
    md_128  data = md_mm_mul_ps(md_mm_cvtepi32_ps(coord), md_mm_set1_ps(inv_precision));
    //md_128i mask = md_mm_set_epi32(0, -1, -1, -1);
    //md_mm_maskstore_ps(dst, mask, data);
    MEMCPY(dst, &data, 3 * sizeof(float));
}

static inline void init_unpack_data_bits(unpack_data_t* data, uint32_t num_of_bits) {
    uint32_t partbits = num_of_bits & 7;
    data->num_of_bits = (uint8_t)num_of_bits;
    data->big_shift = 64 - ((num_of_bits + 7) & ~7);  // Align to the next multiple of 8
    data->sml_shift = (8 - partbits) & 7;
    data->part_mask = partbits ? (1ul << partbits) - 1 : 0xFF;
}

static inline v4i_t unpack_coord64(br_t* r, const unpack_data_t* unpack) {
    uint64_t w = br_read(r, unpack->num_of_bits);
    uint64_t t = ((w << unpack->sml_shift) & (~0xFFull)) | (w & unpack->part_mask);
    uint64_t v = BSWAP64(t) >> unpack->big_shift;
    uint32_t x = (uint32_t)DIV(v, &unpack->div_zy);
    uint64_t q = v - x * unpack->size_zy;
    uint32_t y = (uint32_t)DIV(q, &unpack->div_z);
    uint32_t z = (uint32_t)(q - y * unpack->size_z);

    return v4i_set(x, y, z, 0);
}

static inline v4i_t unpack_coords64_2(uint64_t data, const unpack_data_t* unpack) {
    uint64_t t = ((data << unpack->sml_shift) & (~0xFFull)) | (data & unpack->part_mask);
    uint64_t v = BSWAP64(t) >> unpack->big_shift;
    uint32_t x = (uint32_t)DIV(v, &unpack->div_zy);
    uint64_t q = v - x * unpack->size_zy;
    uint32_t y = (uint32_t)DIV(q, &unpack->div_z);
    uint32_t z = (uint32_t)(q - y * unpack->size_z);

    return v4i_set(x, y, z, 0);
}

static inline v4i_t unpack_coord128(br_t* r, const unpack_data_t* unpack) {
    int fullbytes = unpack->num_of_bits >> 3;
    int partbits  = unpack->num_of_bits  & 7;
    ASSERT(fullbytes >= 8);

#ifdef HAS_INT128_T
    __uint128_t v = BSWAP64(br_read(r, 64));
    int i = 8;
    for (; i < fullbytes; i++) {
        v |= ((__uint128_t) br_read(r, 8)) << (8 * i);
    }
    if (partbits) {
        v |= ((__uint128_t) br_read(r, partbits)) << (8 * i);
    }
    uint32_t x = (uint32_t)(v / unpack->size_zy);
    uint64_t q = (uint64_t)(v - x*unpack->size_zy);
    uint32_t y = (uint32_t)(q / unpack->size_z);
    uint32_t z = (uint32_t)(q % unpack->size_z);
    return v4i_set(x, y, z, 0);
#elif defined(__x86_64__) && defined(_MSC_VER) && (_MSC_VER >= 1920)
    uint64_t lo = BSWAP64(br_read(r, 64));
    uint64_t hi = 0;
    fullbytes -= 8;
    int i = 0;
    for (; i < fullbytes; ++i) {
        hi |= ((uint64_t)br_read(r, 8)) << (8 * i);
    }
    if (partbits) {
        hi |= ((uint64_t)br_read(r, partbits)) << (8 * i);
    }

    uint64_t q;
    uint32_t x = (uint32_t)_udiv128(hi, lo, unpack->size_zy, &q);
    uint32_t y = (uint32_t)(q / unpack->size_z);
    uint32_t z = (uint32_t)(q % unpack->size_z);
    return v4i_set(x, y, z, 0);
#else
    // Default fallback
    uint32_t nums[4] = {0};
    uint32_t sizes[3] = {0, (uint32_t)(unpack->size_zy / unpack->size_z), unpack->size_z};
    uint32_t bytes[16];
    bytes[1] = bytes[2] = bytes[3] = 0;
    int num_of_bytes = 0;
    for (int i = 0; i < fullbytes; ++i) {
        bytes[num_of_bytes++] = (uint32_t)br_read(r, 8);
    }
    if (partbits) {
        bytes[num_of_bytes++] = (uint32_t)br_read(r, partbits);
    }
    for (int i = 2; i > 0; i--) {
        uint32_t num = 0;
        for (int j = num_of_bytes - 1; j >= 0; j--) {
            num = (num << 8) | bytes[j];
            bytes[j] = num / sizes[i];
            num = num % sizes[i];
        }
        nums[i] = num;
    }
    nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
    return v4i_load(nums);
#endif
}

static inline v4i_t br_read_ints(br_t* br, const uint32_t bitsizes[]) {
    uint32_t val[4];
    for (int i = 0; i < 3; ++i) {
        uint32_t num_of_bits = bitsizes[i];
        uint32_t partbits  = num_of_bits & 7;
        uint32_t big_shift = 32 - ((num_of_bits + 7) & ~7);  // Align to the next multiple of 8
        uint32_t sml_shift = (8 - partbits) & 7;
        uint32_t part_mask = partbits ? (1 << partbits) - 1 : 0xFF;
        uint32_t w = (uint32_t)br_read(br, num_of_bits);
        uint32_t r = ((w << sml_shift) & (~0xFF)) | (w & part_mask);
        uint32_t v = BSWAP32(r) >> big_shift;
        val[i] = v;
    }
    return v4i_load(val);
}

static inline __m256i _mm256_bswap_epi64(__m256i x) {
    const __m256i shuffle_mask = _mm256_set_epi8(
        8, 9,10,11,12,13,14,15,
        0, 1, 2, 3, 4, 5, 6, 7,
        8, 9,10,11,12,13,14,15,
        0, 1, 2, 3, 4, 5, 6, 7
    );
    return _mm256_shuffle_epi8(x, shuffle_mask);
}

static inline __m128i _mm_bswap_epi64(__m128i x) {
    const __m128i shuffle_mask = _mm_set_epi8(
        8, 9,10,11,12,13,14,15,
        0, 1, 2, 3, 4, 5, 6, 7
    );
    return _mm_shuffle_epi8(x, shuffle_mask);
}


static inline __m256i load_and_extract_bits_avx2_4_56_le(const uint8_t* data, int bit_offset, const unpack_data_t* unpack) {
    // Compute bit offsets per lane
    __m256i idx         = md_mm256_set_epi64(3,2,1,0);
    __m256i base_bit    = md_mm256_set1_epi64(bit_offset);
    __m256i bit_offsets = _mm256_add_epi64(base_bit, _mm256_mul_epu32(idx, md_mm256_set1_epi64(unpack->num_of_bits)));

    // Byte offset for each lane
    __m256i byte_offsets = _mm256_srli_epi64(bit_offsets, 3);

    // Bit shifts within each 32-bit load
    __m256i shift       = _mm256_and_si256(bit_offsets, md_mm256_set1_epi64(7));

    // Load 8 bytes from base
    __m256i x           = _mm256_i64gather_epi64((const __int64_t*)data, byte_offsets, 1);

    __m256i y          = _mm256_bswap_epi64(x);

    __m256i hi_shift = _mm256_sub_epi64(_mm256_set1_epi64x(unpack->hi_shift), shift);
    __m256i lo_shift = _mm256_sub_epi64(_mm256_set1_epi64x(unpack->lo_shift), shift);

    __m256i hi_mask  = _mm256_set1_epi64x(unpack->hi_mask);
    __m256i lo_mask  = _mm256_set1_epi64x(unpack->lo_mask);

    __m256i hi = _mm256_and_si256(_mm256_srlv_epi64(y, hi_shift), hi_mask);
    __m256i lo = _mm256_and_si256(_mm256_srlv_epi64(y, lo_shift), lo_mask);

    __m256i combined = _mm256_or_si256(hi, lo);

    __m256i swapped  = _mm256_bswap_epi64(combined);

    __m256i result   = _mm256_srlv_epi64(swapped, md_mm256_set1_epi64(unpack->big_shift));
    return result;
}

static inline __m256i extract_bits_avx2_4_64_le(__m256i data, __m256i shift, const unpack_data_t* unpack) {
    __m256i data_be  = _mm256_bswap_epi64(data);

    __m256i hi_shift = _mm256_sub_epi64(_mm256_set1_epi64x(unpack->hi_shift), shift);
    __m256i lo_shift = _mm256_sub_epi64(_mm256_set1_epi64x(unpack->lo_shift), shift);

    __m256i hi_mask  = _mm256_set1_epi64x(unpack->hi_mask);
    __m256i lo_mask  = _mm256_set1_epi64x(unpack->lo_mask);

    __m256i hi       = _mm256_and_si256(_mm256_srlv_epi64(data_be, hi_shift), hi_mask);
    __m256i lo       = _mm256_and_si256(_mm256_srlv_epi64(data_be, lo_shift), lo_mask);

    __m256i combined = _mm256_or_si256(hi, lo);
    __m256i swapped  = _mm256_bswap_epi64(combined);

    __m256i result   = _mm256_srlv_epi64(swapped, md_mm256_set1_epi64(unpack->big_shift));
    return result;
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
#define USE_STATELESS 1

static inline __m256i create_correction_mask(int8_t cor_a, int8_t cor_b) {
    #if 0
    __m128i broad_a = _mm_set1_epi8(cor_a);
    __m128i broad_b = _mm_set1_epi8(cor_b);
    __m256i combined = _mm256_set_m128i(broad_b, broad_a);
    __m256i mask = _mm256_setr_epi8(
        0,0,0,0,0,0,0,0,
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF, // 0
        0,0,0,0,0,0,0,0,
        0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF);  // 0
    return _mm256_and_si256(combined, mask);
    #else
    __m256i broad_a = _mm256_set1_epi8(cor_a);
    __m256i broad_b = _mm256_set1_epi8(cor_b);

    __m256i mask = _mm256_setr_epi8(
        // lower 128-bit lane
        0, 0, 0, 0, 0, 0, 0, 0,    // use broad_b
        (char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80, // zero
        // upper 128-bit lane
        0, 0, 0, 0, 0, 0, 0, 0,    // use broad_a
        (char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80  // zero
    );

    // blend broad_a and broad_b using shuffle mask
    __m256i blended = _mm256_blendv_epi8(broad_b, broad_a, _mm256_setr_epi8(
        // lower 128-bit
        0,0,0,0,0,0,0,0,    // b
        (char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF, // 0
        // upper 128-bit
        (char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF, // a
        (char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF  // 0
    ));

    // zero out the zeros
    __m256i result = _mm256_shuffle_epi8(blended, mask);
    return result;
    #endif
}

size_t md_xtc_read_frame_coords(md_file_o* xdr_file, float* out_coords_ptr, size_t out_coord_cap) {
    if (xdr_file == NULL || out_coords_ptr == NULL) {
        return 0;
    }

    int natoms;
    if (!xdr_read_int32(&natoms, 1, xdr_file)) {
        return 0;
    }

    if (out_coord_cap < (size_t)natoms) {
        MD_LOG_ERROR("The supplied coordinate buffer is not sufficient to fit the number of atoms within the frame");
        return 0;
    }

    /* Dont bother with compression for three atoms or less */
    if (natoms <= 9) {
        if (!xdr_read_float(out_coords_ptr, natoms * 3, xdr_file)) {
            return 0;
        }
        return (size_t)natoms;
    }

    /* Compression-time if we got here. Read precision first */
    float precision;
    if (!xdr_read_float(&precision, 1, xdr_file)) {
        return 0;
    }

    int32_t minint[3], maxint[3];
    if (!xdr_read_int32(minint, 3, xdr_file) || !xdr_read_int32(maxint, 3, xdr_file)) {
        return 0;
    }

    uint32_t sizeint[3] = {
        (uint32_t)(maxint[0] - minint[0] + 1),
        (uint32_t)(maxint[1] - minint[1] + 1),
        (uint32_t)(maxint[2] - minint[2] + 1),
    };

    uint32_t bitsize = 0;
    uint32_t bitsizeint[3] = {0};
    unpack_data_t big_unpack = {0};
    uint32_t big_skip = 0;

    /* check if one of the sizes is to big to be multiplied */
    if ((sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffff) {
        bitsizeint[0] = sizeofint(sizeint[0]);
        bitsizeint[1] = sizeofint(sizeint[1]);
        bitsizeint[2] = sizeofint(sizeint[2]);
        big_skip = bitsizeint[0] + bitsizeint[1] + bitsizeint[2];
    } else {
        bitsize = sizeofints(3, sizeint);
        big_unpack.size_zy = (uint64_t)sizeint[1] * sizeint[2];
        big_unpack.size_z = sizeint[2];
        big_unpack.num_of_bits = (uint8_t)bitsize;
        big_unpack.div_zy = DIV_INIT(big_unpack.size_zy);
        big_unpack.div_z = DIV_INIT(big_unpack.size_z);
        init_unpack_data_bits(&big_unpack, bitsize);
        big_skip = bitsize;
    }

    int smallidx;
    if (!xdr_read_int32(&smallidx, 1, xdr_file)) {
        return 0;
    }

    int idx = MAX(smallidx - 1, FIRSTIDX);
    int smaller = magicints[idx] / 2;
    int smallnum = magicints[smallidx] / 2;
    uint32_t smallsize = magicints[smallidx];

    unpack_data_t sml_unpack = {
        .size_zy = (uint64_t)smallsize * smallsize,
        .size_z  = smallsize,
        .num_of_bits = (uint8_t)smallidx,
        .div_zy = denoms_64_2[smallidx - FIRSTIDX],
        .div_z  = denoms_64_1[smallidx - FIRSTIDX],
    };
    init_unpack_data_bits(&sml_unpack, smallidx);

    /* length in bytes */
    uint32_t num_bytes = 0;
    if (!xdr_read_int32((int32_t*)&num_bytes, 1, xdr_file)) {
        return 0;
    }

    size_t temp_pos = md_temp_get_pos();
    size_t mem_size = ALIGN_TO(num_bytes, 32);
    uint8_t* mem = md_temp_push_aligned(mem_size, 32);

    int atom_idx = 0;
    if (!xdr_read_bytes(mem, num_bytes, xdr_file)) {
        goto done;
    }

    const uint8_t* base = mem;
    size_t bit_offset = 0;

    br_t br = {0};
    br_init(&br, mem, mem_size / 8);

    float* lfp = out_coords_ptr;
    float inv_precision = 1.0f / precision;
    v4i_t vminint = v4i_set(minint[0], minint[1], minint[2], 0);
    v4i_t thiscoord;
    int run = 0;
    int run_count = 0;

#ifdef PRINT_DEBUG
    uint32_t smlbits_hist[LASTIDX] = {0};
    uint32_t rle_hist[9] = {0};
    MD_LOG_DEBUG("bigsize: %i", bitsize);
#endif

    size_t bit_offset_lo = big_skip;
    size_t bit_offset_hi = big_skip + 4 * sml_unpack.num_of_bits;

    while (atom_idx < natoms) {

#if 1
        size_t offset = bit_offset >> 3;
        size_t shift  = bit_offset  & 7;
        __m128i v = _mm_loadu_si128(base + offset);
        uint64_t x = _mm_cvtsi128_si64(v);

        x = BSWAP64(x);
        x <<= shift;

        uint64_t coord_bits = x >> (64 - big_unpack.num_of_bits);
        uint64_t data       = x >> (64 - big_unpack.num_of_bits - 6) & 63;
        uint32_t flag       = data & 32;
        uint32_t skip       = flag ? 6 : 1;
        uint32_t sml_shift  = flag ? 5 : 0;
        thiscoord = unpack_coords64_2(coord_bits, &big_unpack);

        size_t bit_off_base = bit_offset + big_skip + 1;

        __m256i idx_lo    = _mm256_set_epi64x(3,2,1,0);
        __m256i idx_hi    = _mm256_set_epi64x(7,6,5,4);
        __m256i base_lo   = _mm256_add_epi64(_mm256_set1_epi64x(bit_off_base), _mm256_mul_epu32(idx_lo, _mm256_set1_epi64x(sml_unpack.num_of_bits)));
        __m256i base_hi   = _mm256_add_epi64(_mm256_set1_epi64x(bit_off_base), _mm256_mul_epu32(idx_hi, _mm256_set1_epi64x(sml_unpack.num_of_bits)));
        __m256i shift_lo  = _mm256_and_si256(base_lo, _mm256_set1_epi64x(7));
        __m256i shift_hi  = _mm256_and_si256(base_hi, _mm256_set1_epi64x(7));

        size_t off_a =  (bit_off_base + 0 * sml_unpack.num_of_bits) >> 3;
        int8_t cor_a = ((bit_off_base + 1 * sml_unpack.num_of_bits) >> 3) - (off_a + 16);
        size_t off_b =  (bit_off_base + 2 * sml_unpack.num_of_bits) >> 3;
        int8_t cor_b = ((bit_off_base + 3 * sml_unpack.num_of_bits) >> 3) - (off_b + 16);
        size_t off_c =  (bit_off_base + 4 * sml_unpack.num_of_bits) >> 3;
        int8_t cor_c = ((bit_off_base + 5 * sml_unpack.num_of_bits) >> 3) - (off_c + 16);
        size_t off_d =  (bit_off_base + 6 * sml_unpack.num_of_bits) >> 3;
        int8_t cor_d = ((bit_off_base + 7 * sml_unpack.num_of_bits) >> 3) - (off_d + 16);

        __m256i data_lo = _mm256_loadu2_m128((const __m128i*)(base + off_b), (const __m128i*)(base + off_a));
        __m256i data_hi = _mm256_loadu2_m128((const __m128i*)(base + off_d), (const __m128i*)(base + off_c));

        //__m256i lo = load_and_extract_bits_avx2_4_56_le(base, offset_lo, &sml_unpack);
        //__m256i hi = load_and_extract_bits_avx2_4_56_le(base, offset_hi, &sml_unpack);

        bit_offset += big_skip + skip;



#else
        if (bitsize == 0) {
            thiscoord = br_read_ints(&br, bitsizeint);
        } else {
            if (big_unpack.num_of_bits <= 64) {
#if USE_STATELESS
                thiscoord = unpack_coords64_2(extract_u64_be(base, bit_offset, big_unpack.num_of_bits), &big_unpack);
#else
                thiscoord = unpack_coord64(&br, &big_unpack);
#endif
            } else {
                thiscoord = unpack_coord128(&br, &big_unpack);
            }
        }
#endif
        thiscoord = v4i_add(thiscoord, vminint);

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

#ifdef PRINT_DEBUG
        rle_hist[run_count]++;
#endif

        if (run > 0) {
            bit_offset += run_count * sml_unpack.num_of_bits;

            __m256i cor_lo = create_correction_mask(cor_a, cor_b);
            __m256i cor_hi = create_correction_mask(cor_c, cor_d);

            shift_lo = _mm256_add_epi64(shift_lo, _mm256_set1_epi64x(sml_shift));
            shift_hi = _mm256_add_epi64(shift_hi, _mm256_set1_epi64x(sml_shift));

            // Construct a shuffle mask to emulate a 64-bit gather on the
            __m256i shuffle_base = _mm256_set_epi8(15,14,13,12,11,10,9,8, 7,6,5,4,3,2,1,0, 15,14,13,12,11,10,9,8, 7,6,5,4,3,2,1,0);
            __m256i shuffle_lo   = _mm256_sub_epi8(shuffle_base, cor_lo);
            __m256i shuffle_hi   = _mm256_sub_epi8(shuffle_base, cor_hi);

            // Extract the correct data from lo and hi
            data_lo = _mm256_shuffle_epi8(data_lo, shuffle_lo);
            data_hi = _mm256_shuffle_epi8(data_hi, shuffle_hi);

            // Create mask for store
            __m256i v_run  = _mm256_set1_epi32(run);
            __m256i mask_a = _mm256_cmpgt_epi32(v_run, _mm256_set_epi32( 7, 6, 5, 4, 3, 2, 1, 0));
            __m256i mask_b = _mm256_cmpgt_epi32(v_run, _mm256_set_epi32(15,14,13,12,11,10, 9, 8));
            __m256i mask_c = _mm256_cmpgt_epi32(v_run, _mm256_set_epi32(23,22,21,20,19,18,17,16));

            //__m256i lo = load_and_extract_bits_avx2_4_56_le(base, bit_offset, &sml_unpack);
            //__m256i hi = load_and_extract_bits_avx2_4_56_le(base, bit_offset + 4 * sml_unpack.num_of_bits, &sml_unpack);
            //lo = extract_bits_avx2_4_64_le(lo, &sml_unpack);
            //hi = extract_bits_avx2_4_64_le(hi, &sml_unpack);

            __m256 x, y, z;
            x = _mm256_mul_ps(_mm256_cvtepi32_ps(data_lo), _mm256_set1_ps(inv_precision));
            y = _mm256_mul_ps(_mm256_cvtepi32_ps(data_lo), _mm256_set1_ps(inv_precision));
            z = _mm256_mul_ps(_mm256_cvtepi32_ps(data_hi), _mm256_set1_ps(inv_precision));

            // Interleave X, Y and Z into 3x XYZXYZ...
            //__m256 a, b, c;

            _mm256_maskstore_ps(lfp +  0, mask_a, x);
            _mm256_maskstore_ps(lfp +  8, mask_b, y);
            _mm256_maskstore_ps(lfp + 16, mask_c, z);

            lfp += run;
        } else {
            write_coord(lfp, thiscoord, inv_precision); lfp += 3;
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
        if (smallidx != sml_unpack.num_of_bits) {
#ifdef PRINT_DEBUG
            smlbits_hist[smallidx]++;
#endif
            uint32_t sml_size       = magicints[smallidx];
            sml_unpack.size_zy      = (uint64_t)sml_size * sml_size;
            sml_unpack.size_z       = sml_size;
            sml_unpack.num_of_bits  = (uint8_t)smallidx;
            sml_unpack.div_zy       = denoms_64_2[smallidx - FIRSTIDX];
            sml_unpack.div_z        = denoms_64_1[smallidx - FIRSTIDX];
            init_unpack_data_bits(&sml_unpack, smallidx);
        }
    }
    
done:
    md_temp_set_pos_back(temp_pos);

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
            header->unit_cell = md_util_unit_cell_from_matrix(box);
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

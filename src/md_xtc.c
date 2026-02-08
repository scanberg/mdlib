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
// --- Normalized sequential bit reader --------------------------------------
// * Input bitstream is big-endian (network order), MSB-first.
// * Returned integers are canonical: bitfield is placed in the low bits of the return value.

typedef struct brn_t {
    const uint64_t* p;      // Next qword to load
    uint64_t a;             // Current qword (BE->host swapped)
    uint64_t b;             // Next qword (BE->host swapped)
    uint32_t bitpos;        // 0..63, number of bits consumed from MSB side of a
    uint32_t qwords_left;   // Remaining qwords in p
} brn_t;

static inline uint64_t brn_load_qword(const uint64_t* src) {
    uint64_t v = *src;
#if __LITTLE_ENDIAN__
    v = BSWAP64(v);
#endif
    return v;
}

static inline void brn_init(brn_t* r, const void* data, size_t num_bytes) {
    ASSERT(r);
    ASSERT(data);
    ASSERT((num_bytes & 7) == 0);
    const uint64_t* q = (const uint64_t*)data;
    const size_t num_qwords = num_bytes / 8;
    ASSERT(num_qwords >= 2);

    r->a = brn_load_qword(q + 0);
    r->b = brn_load_qword(q + 1);
    r->p = q + 2;
    r->qwords_left = (uint32_t)num_qwords - 2;
    r->bitpos = 0;
}

static inline void brn_slide_1(brn_t* r) {
    r->a = r->b;
    if (r->qwords_left) {
        r->b = brn_load_qword(r->p++);
        r->qwords_left -= 1;
    } else {
        r->b = 0;
    }
}

// Produce a 64-bit "view" of the stream starting at current bitpos (MSB-first).
static inline uint64_t brn_view64(const brn_t* r) {
    const uint32_t s = r->bitpos;

    // view = (a << s) | (b >> (64-s)) but avoid UB when s==0.
    const uint64_t upper = r->a << s;
    const uint64_t nz_mask = -(uint64_t)(s != 0);
    const uint64_t lower = (r->b >> (64 - s)) & nz_mask;
    return upper | lower;
}

static inline uint64_t brn_peek_u64(const brn_t* r, uint32_t n) {
    ASSERT(r);
    ASSERT(n >= 1 && n <= 64);
    const uint64_t view = brn_view64(r);
    return view >> (64 - n);
}

static inline void brn_skip(brn_t* r, uint32_t n) {
    ASSERT(r);
    ASSERT(n >= 1 && n <= 64);
    const uint32_t newpos = r->bitpos + n;
    r->bitpos = newpos & 63u;
    if (newpos >= 64) {
        brn_slide_1(r);
        if (newpos >= 128) {
            brn_slide_1(r);
        }
    }
}

// Read 1..64 bits, return in low bits.
static inline uint64_t brn_read_u64(brn_t* r, uint32_t n) {
    ASSERT(r);
    ASSERT(n >= 1 && n <= 64);

    const uint64_t view = brn_view64(r);
    const uint64_t res  = view >> (64 - n);

    const uint32_t newpos = r->bitpos + n;
    r->bitpos = newpos & 63u;

    if (newpos >= 64) {
        brn_slide_1(r);
        // n can be 64 and bitpos can be non-zero: could cross 2 qwords.
        if (newpos >= 128) {
            brn_slide_1(r);
        }
    }

    return res;
}

// Read 65..128 bits into canonical (hi, lo) where V = (hi<<64)|lo and uses low n bits.
static inline void brn_read_u128(brn_t* r, uint32_t n, uint64_t out_hi_lo[2]) {
    ASSERT(r);
    ASSERT(out_hi_lo);
    ASSERT(n > 64 && n <= 128);

    // Strategy:
    // - Read 64 bits first (high part of the n-bit field)
    // - Read remaining (n-64) bits (low part)
    //
    // But we must preserve bit-exact order. The field is MSB-first.
    //
    // So: hi_part contains the first 64 bits of the field.
    //      lo_part contains the remaining bits left-aligned within 64, then shift it down.
    const uint32_t lo_bits = n - 64;

    const uint64_t hi_part = brn_read_u64(r, 64);
    const uint64_t lo_part = brn_read_u64(r, lo_bits);

    // Now normalize: field should occupy the low n bits of a 128-bit integer.
    // hi_part currently holds 64 bits, lo_part holds lo_bits bits.
    // Place them as: V = (hi_part << lo_bits) | lo_part.
    // That means:
    //   out_hi = hi_part >> (64 - lo_bits) ??? No, because hi_part is 64 bits already.
    //
    // Compute in 128-bit-free way:
    if (lo_bits == 64) {
        out_hi_lo[0] = hi_part;
        out_hi_lo[1] = lo_part;
    } else {
        // V = ( (__uint128)hi_part << lo_bits ) | lo_part
        out_hi_lo[0] = (lo_bits == 0) ? 0 : (hi_part >> (64 - lo_bits));
        out_hi_lo[1] = (hi_part << lo_bits) | lo_part;
    }
}

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

static inline md_128 cvt_coord(v4i_t v, md_128 invp) {
	return md_mm_mul_ps(md_mm_cvtepi32_ps(v), invp);
}

static inline void write_coord_interleaved(float* dst, md_128 coord) {
    MEMCPY(dst, &coord, sizeof(float) * 3);
}

static inline void write_coord_xyz(float* x, float* y, float* z, int offset, md_128 coord) {
    float tmp[4];
    md_mm_storeu_ps(tmp, coord);
    x[offset] = tmp[0];
    y[offset] = tmp[1];
    z[offset] = tmp[2];
}

static inline void init_unpack_bit_data(bit_data_t* data, uint32_t num_of_bits) {
    uint32_t partbits = num_of_bits & 7;
    data->num_of_bits = (uint8_t)num_of_bits;
    data->big_shift = 64 - ((num_of_bits + 7) & ~7);  // Align to the next multiple of 8
    data->sml_shift = (8 - partbits) & 7;
    data->part_mask = partbits ? (1ul << partbits) - 1 : 0xFF;
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

    uint64_t raw[2];
    MEMCPY(raw, base + byte_offset, sizeof(raw));

#if __LITTLE_ENDIAN__
    raw[0] = BSWAP64(raw[0]);
    raw[1] = BSWAP64(raw[1]);
#endif

    uint64_t hi = (raw[0] << bit_in_byte);
    uint64_t lo = (raw[1] << bit_in_byte) | (raw[0] >> (64 - bit_in_byte));
    out[0] = hi;
    out[1] = lo >> shift_right;
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
    
#if __LITTLE_ENDIAN__
    raw[0] = BSWAP64(raw[0]);
    raw[1] = BSWAP64(raw[1]);
#endif

	uint64_t hi = raw[0];
	uint64_t lo = raw[1];
    
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

#if __LITTLE_ENDIAN__
    raw = BSWAP64(raw);
#endif

    // Extract the bits (they're now at the top after bswap)
	// Mask upper by shifting left and then right
    return (raw << bit_in_byte) >> (64 - bit_length);
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
#if __LITTLE_ENDIAN__
    raw = BSWAP64(raw);
#endif
    return (uint32_t)((raw << bit_in_byte) >> (64 - bit_length));
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
#if __LITTLE_ENDIAN__
    raw = BSWAP32(raw);
#endif
    return (uint32_t)((raw << bit_in_byte) >> (32 - bit_length));
}

// Generic version that issues the appropriate version based on number of bits
static inline void extract_bits_be(
    uint64_t out[2],
    const uint8_t* base,
    size_t bit_offset,
    size_t bit_length)   // 1..121
{
    if (bit_length <= 57) {
        out[0] = extract_bits_be_raw_57(base, bit_offset, bit_length);
    } else if (bit_length <= 64) {
        out[0] = extract_bits_be_raw_64(base, bit_offset, bit_length);
    } else {
        extract_bits_be_raw_121(out, base, bit_offset, bit_length);
    }
}

static inline v4i_t unpack_coord64(uint64_t w, const unpack_data_t* unpack) {
	// Reconstruct correct value with shifts and masks
    uint64_t t = ((w << unpack->bit.sml_shift) & (~0xFFull)) | (w & unpack->bit.part_mask);
    uint64_t v = BSWAP64(t) >> unpack->bit.big_shift;

    uint32_t x  = (uint32_t)DIV(v, &unpack->div_zy);   // v / (Y*Z)
    uint64_t yz = (uint64_t)DIV(v, &unpack->div_z);    // v / Z

    uint32_t y = (uint32_t)(yz - (uint64_t)x  * unpack->size_y);
    uint32_t z = (uint32_t)(v  - (uint64_t)yz * unpack->size_z);

	return v4i_set(x, y, z, 0);
}

static inline v4i_t unpack_coord128(uint64_t w[2], const unpack_data_t* unpack) {
    const uint64_t zy = (uint64_t)unpack->size_z * unpack->size_y;

	// Construct correct hi and lo portions of the 128-bit integer
	// One half is just a BSWAP, the other needs shifting and masking as seen in unpack_coord64
    uint64_t hi = ((w[1] << unpack->bit.sml_shift) & (~0xFFull)) | (w[1] & unpack->bit.part_mask);
	hi = BSWAP64(hi) >> unpack->bit.big_shift;
	uint64_t lo = BSWAP64(w[0]);

#ifdef HAS_INT128_T
	__uint128_t v = (__uint128_t)lo | (__uint128_t)hi << 64;
    uint32_t x = (uint32_t)(v / zy);
    uint64_t q = (uint64_t)(v - x*zy);
    uint32_t y = (uint32_t)(q / unpack->size_z);
    uint32_t z = (uint32_t)(q % unpack->size_z);
    return v4i_set(x, y, z, 0);
#elif defined(__x86_64__) && defined(_MSC_VER) && (_MSC_VER >= 1920)
    uint64_t q = 0;
    uint32_t x = (uint32_t)_udiv128(hi, lo, zy, &q);
    uint32_t y = (uint32_t)(q / unpack->size_z);
    uint32_t z = (uint32_t)(q % unpack->size_z);
    return v4i_set(x, y, z, 0);

#else
    // Default fallback
    // limbs[0] is least-significant 32 bits

    const int fullbytes = unpack->bit.num_of_bits >> 3;
    const int partbits  = unpack->bit.num_of_bits  & 7;
    const int num_of_bytes = fullbytes + (partbits ? 1 : 0);
    const uint32_t sizes[3] = {0, unpack->size_y, unpack->size_z};
    uint32_t nums[4]  = {0};
    uint32_t limbs[4] = {0};
    
	// Load limbs
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

// Loads an unaligned BE qword and returns it as a host uint64_t with BE->host byte swap.
static inline uint64_t load_be_u64_unaligned(const uint8_t* p) {
    uint64_t v;
    MEMCPY(&v, p, sizeof(v));
#if __LITTLE_ENDIAN__
    v = BSWAP64(v);
#endif
    return v;
}

// Extract 1..64 bits starting at bit_offset (MSB-first), return normalized in low bits.
// ASSUMES input stream is big-endian bit/byte order (as your existing extractors do).
static inline uint64_t extract_bits_be_norm_64(const uint8_t* base, size_t bit_offset, uint32_t bit_length) {
    ASSERT(base);
    ASSERT(bit_length >= 1 && bit_length <= 64);

    const size_t byte_offset = bit_offset >> 3;
    const uint32_t bit_in_byte = (uint32_t)(bit_offset & 7);

    // We align to a BE 64-bit chunk starting at byte_offset.
    // If the requested bits fit within this qword after applying bit_in_byte, we only need 8 bytes.
    // Condition: bit_in_byte + bit_length <= 64
    if (bit_in_byte + bit_length <= 64) {
        const uint64_t hi = load_be_u64_unaligned(base + byte_offset);

        // Treat hi as MSB-first shift register: shift left by bit_in_byte then right to low bits.
        return (hi << bit_in_byte) >> (64 - bit_length);
    } else {
        // Need two qwords
        const uint64_t hi = load_be_u64_unaligned(base + byte_offset);
        const uint64_t lo = load_be_u64_unaligned(base + byte_offset + 8);

        // Build a 64-bit view starting at bit_in_byte across the 128-bit window (hi||lo), MSB-first.
        const uint64_t view = (hi << bit_in_byte) | (lo >> (64 - bit_in_byte));

        return view >> (64 - bit_length);
    }
}

static inline v4i_t extract_and_unpack(
    const uint8_t* base,
    size_t bit_offset,
    size_t bit_length,
    const unpack_data_t* unpack)
{
    if (bit_length <= 64) {
        uint64_t w;
        if (bit_length <= 57)
			w = extract_bits_be_raw_57(base, bit_offset, bit_length);
        else
            w = extract_bits_be_raw_64(base, bit_offset, bit_length);
        return unpack_coord64(w, unpack);
    } else {
        uint64_t w[2];
        extract_bits_be_raw_121(w, base, bit_offset, bit_length);
        return unpack_coord128(w, unpack);
	}
}

static inline uint32_t unpack_uint32(uint32_t w, const bit_data_t* bits) {
	uint32_t t = ((w << bits->sml_shift) & (~0xFFul)) | (w & bits->part_mask);
	uint32_t v = BSWAP32(t) >> bits->big_shift;
    return v;
}

static inline v4i_t extract_ints3(const uint8_t* base, size_t bit_offset, const bit_data_t bits[3]) {
    uint32_t v[4];
    for (int i = 0; i < 3; ++i) {
        size_t num_bits = bits[i].num_of_bits;
        uint32_t w = extract_bits_be_raw_32(base, bit_offset, num_bits);
        v[i] = unpack_uint32(w, &bits[0]);
        bit_offset += num_bits;
    }
    return v4i_load(v);
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
#define USE_BITREADER 1

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
    if (!xdr_read_int32(&smallidx, 1, xdr_file)) {
        return 0;
    }

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
    if (!xdr_read_int32((int32_t*)&num_bytes, 1, xdr_file)) {
        return 0;
    }
    size_t mem_size = ALIGN_TO(num_bytes, 32);
    void* mem = malloc(mem_size);

    uint8_t* base = mem;
    size_t bit_offset = 0;

    int atom_idx = 0;
    if (!xdr_read_bytes((uint8_t*)mem, num_bytes, xdr_file)) {
        goto done;
    }

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

    while (atom_idx < natoms) {
        uint32_t run_data = extract_bits_be_raw_25(base, bit_offset + big_bits, 6);

        if (bitsize == 0) {
            thiscoord = extract_ints3(base, bit_offset, big_bit_data);
        } else {
			thiscoord = extract_and_unpack(base, bit_offset, bitsize, &big_unpack);
        }
        thiscoord = v4i_add(thiscoord, vminint);

        uint32_t run_flag = run_data & 32;
        uint32_t run_bits = run_flag ? 6 : 1;

        bit_offset += big_bits + run_bits;

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

#ifdef PRINT_DEBUG
        rle_hist[run_count]++;
#endif

        // If run, write first atom after second atom
        int offset = run ? 3 : 0;
        write_coord_interleaved(lfp + offset, cvt_coord(thiscoord, invp)); lfp += 3;

        if (run) {
            const v4i_t vsmall = v4i_set1(smallnum);
            const int num_bits = sml_unpack.bit.num_of_bits;

			v4i_t delta = v4i_sub(thiscoord, vsmall);
			v4i_t coord = extract_and_unpack(base, bit_offset, num_bits, &sml_unpack);
            thiscoord = v4i_add(coord, delta);

            // Write second atom before first atom
            write_coord_interleaved(lfp - 3, cvt_coord(thiscoord, invp)); lfp += 3;
            bit_offset += num_bits;

            for (int i = 1; i < run_count; i += 1) {
				delta = v4i_sub(thiscoord, vsmall);
				coord = extract_and_unpack(base, bit_offset, num_bits, &sml_unpack);
                thiscoord = v4i_add(coord, delta);
                write_coord_interleaved(lfp, cvt_coord(thiscoord, invp)); lfp += 3;
                bit_offset += num_bits;
            }
            #ifdef PRINT_DEBUG
                smlbits_hist[smallidx]++;
            #endif
        }

        if (is_smaller) {
            smallidx += is_smaller;
            smallnum = (int)magicints[smallidx] / 2;
            smaller  = (smallidx > FIRSTIDX) ? (int)magicints[smallidx - 1] / 2 : 0;
        }

        if (smallidx != sml_unpack.bit.num_of_bits) {
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

        atom_idx += batch_size;
    }
    
done:
    free(mem);

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

size_t md_xtc_read_frame_coords_xyz(md_file_o* xdr_file, float* out_x, float* out_y, float* out_z, size_t out_coord_cap) {
    if (xdr_file == NULL || out_x == NULL || out_y == NULL || out_z == NULL) {
		MD_LOG_ERROR("Invalid arguments: xdr_file and output coordinate buffers must not be NULL");
        return 0;
    }

    int natoms;
    if (!xdr_read_int32(&natoms, 1, xdr_file)) {
		MD_LOG_ERROR("XTC: Failed to read number of atoms");
        return 0;
    }

    if (out_coord_cap < (size_t)natoms) {
        MD_LOG_ERROR("The supplied coordinate buffer is not sufficient to fit the number of atoms within the frame");
        return 0;
    }

    /* Dont bother with compression for three atoms or less */
    if (natoms <= 9) {
		float coords[9 * 3];
        if (!xdr_read_float(coords, natoms * 3, xdr_file)) {
			MD_LOG_ERROR("XTC: Failed to read uncompressed coordinates");
            return 0;
        }
		// De-interleave coordinates into separate x, y, z arrays
        for (int i = 0; i < natoms; ++i) {
            out_x[i] = coords[i * 3 + 0];
            out_y[i] = coords[i * 3 + 1];
            out_z[i] = coords[i * 3 + 2];
		}
        return (size_t)natoms;
    }

    // Compression-time if we got here. Read precision first
    float precision;
    if (!xdr_read_float(&precision, 1, xdr_file)) {
		MD_LOG_ERROR("XTC: Failed to read precision");
        return 0;
    }

    int32_t minint[3], maxint[3];
    if (!xdr_read_int32(minint, 3, xdr_file) || !xdr_read_int32(maxint, 3, xdr_file)) {
        MD_LOG_ERROR("XTC: Failed to read min/max coordinates");
        return 0;
    }

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
    if (!xdr_read_int32(&smallidx, 1, xdr_file)) {
        return 0;
    }

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
    if (!xdr_read_int32((int32_t*)&num_bytes, 1, xdr_file)) {
        return 0;
    }
    size_t mem_size = ALIGN_TO(num_bytes, 32);
    void* mem = malloc(mem_size);

    uint8_t* base = mem;
    size_t bit_offset = 0;

    int atom_idx = 0;
    if (!xdr_read_bytes((uint8_t*)mem, num_bytes, xdr_file)) {
        goto done;
    }

	const md_128 invp = md_mm_set1_ps(1.0f / precision);
    const v4i_t vminint = v4i_set(minint[0], minint[1], minint[2], 0);
    int run = 0;
    int run_count = 0;

    while (atom_idx < natoms) {
        uint32_t run_data = extract_bits_be_raw_25(base, bit_offset + big_bits, 6);

        v4i_t thiscoord;
        if (bitsize == 0) {
            thiscoord = extract_ints3(base, bit_offset, big_bit_data);
        } else {
			thiscoord = extract_and_unpack(base, bit_offset, bitsize, &big_unpack);
        }
        thiscoord = v4i_add(thiscoord, vminint);

        uint32_t run_flag = run_data & 32;
        uint32_t run_bits = run_flag ? 6 : 1;

        bit_offset += big_bits + run_bits;

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

        const int base_idx = atom_idx;
        if (run) {
            const v4i_t vsmall = v4i_set1(smallnum);
            const int num_bits = sml_unpack.bit.num_of_bits;

			// If run, write first atom after second atom
            write_coord_xyz(out_x, out_y, out_z, base_idx + 1, cvt_coord(thiscoord, invp));

			v4i_t delta = v4i_sub(thiscoord, vsmall);
			v4i_t coord = extract_and_unpack(base, bit_offset, num_bits, &sml_unpack);
            thiscoord = v4i_add(coord, delta);

            // Write second atom before first atom
            write_coord_xyz(out_x, out_y, out_z, base_idx, cvt_coord(thiscoord, invp));
            bit_offset += num_bits;

			for (int i = 2; i < batch_size; i += 1) {
				delta = v4i_sub(thiscoord, vsmall);
				coord = extract_and_unpack(base, bit_offset, num_bits, &sml_unpack);
                thiscoord = v4i_add(coord, delta);
                write_coord_xyz(out_x, out_y, out_z, base_idx + i, cvt_coord(thiscoord, invp));
                bit_offset += num_bits;
            }
            atom_idx = base_idx + batch_size;
        } else {
			write_coord_xyz(out_x, out_y, out_z, base_idx, cvt_coord(thiscoord, invp));
            atom_idx = base_idx + 1;
        }

        if (is_smaller) {
            smallidx += is_smaller;
            smallnum = (int)magicints[smallidx] / 2;
            smaller  = (smallidx > FIRSTIDX) ? (int)magicints[smallidx - 1] / 2 : 0;
        }

        if (smallidx != sml_unpack.bit.num_of_bits) {
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
    }
    
done:
    free(mem);
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
            size_t written_count = md_xtc_read_frame_coords_xyz(file, x, y, z, natoms);
            result = (written_count == natoms);
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

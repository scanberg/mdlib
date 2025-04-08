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

#include <core/md_lru_cache.inl>

#include <xdrfile.h>

#define LIBDIVIDE_AVX2
#include <libdivide.h>

#include <string.h>
#include <stdio.h>

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
    XDRFILE* file;
    md_array(int64_t) frame_offsets;
    md_trajectory_header_t header;
    md_allocator_i* allocator;
    md_mutex_t mutex;    
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

static inline bool xdr_read_int16(int16_t* ptr, size_t count, md_file_o* xdr_file) {
    ASSERT(ptr);
    ASSERT(xdr_file);

    size_t bytes = count * sizeof(int16_t);
    if (md_file_read(xdr_file, ptr, bytes) == bytes) {
#if __LITTLE_ENDIAN__
        for (size_t i = 0; i < count; ++i) {
            ptr[i] = BSWAP16(ptr[i]);
        }
#endif
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

    while (size >= num && num_of_bits < 32) {
        num_of_bits++;
        num <<= 1;
    }
    return num_of_bits;
}

static inline int sizeofints(int num_of_ints, unsigned int sizes[]) {
    int i, num;
    unsigned int num_of_bytes, num_of_bits, bytes[32], bytecnt, tmp;
    num_of_bytes = 1;
    bytes[0] = 1;
    num_of_bits = 0;
    for (i = 0; i < num_of_ints; i++) {
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

static inline int decodebits(int buf[3], int num_of_bits) {
    int cnt, num;
    unsigned int lastbits, lastbyte;
    unsigned char* cbuf;
    int mask = (1 << num_of_bits) - 1;

    cbuf = ((unsigned char*)buf) + 3 * sizeof(*buf);
    cnt = buf[0];
    lastbits = (unsigned int)buf[1];
    lastbyte = (unsigned int)buf[2];

    num = 0;
    while (num_of_bits >= 8) {
        lastbyte = (lastbyte << 8) | cbuf[cnt++];
        num |= (lastbyte >> lastbits) << (num_of_bits - 8);
        num_of_bits -= 8;
    }
    if (num_of_bits > 0) {
        if (lastbits < num_of_bits) {
            lastbits += 8;
            lastbyte = (lastbyte << 8) | cbuf[cnt++];
        }
        lastbits -= num_of_bits;
        num |= (lastbyte >> lastbits) & ((1 << num_of_bits) - 1);
    }
    num &= mask;
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    return num;
}

/*
* decodeints - decode 'small' integers from the buf array
*
* this routine is the inverse from encodeints() and decodes the small integers
* written to buf by calculating the remainder and doing divisions with
* the given sizes[]. You need to specify the total number of bits to be
* used from buf in num_of_bits.
*
*/

static inline void decodeints(int buf[], int num_of_ints, int num_of_bits, unsigned int sizes[], int nums[3]) {

    int bytes[16];
    int i, j, num_of_bytes, p, num, size;

    bytes[1] = bytes[2] = bytes[3] = 0;
    num_of_bytes = 0;
    while (num_of_bits > 8) {
        bytes[num_of_bytes++] = decodebits(buf, 8);
        num_of_bits -= 8;
    }
    if (num_of_bits > 0) {
        bytes[num_of_bytes++] = decodebits(buf, num_of_bits);
    }
    //printf("num_bytes %i\n", num_of_bytes);
    for (i = num_of_ints - 1; i > 0; i--) {
        num = 0;
        size = sizes[i];
        for (j = num_of_bytes - 1; j >= 0; j--) {
            num = (num << 8) | bytes[j];
            p = num / size;
            num = num - p * size;
            bytes[j] = p;
        }
        nums[i] = num;
    }
    nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
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

static const struct libdivide_u32_branchfree_t denoms32[] = {
    {0, 2}, {0X9999999A, 3}, {0X55555556, 3}, {0, 3}, {0X9999999A, 4}, {0X47AE147B, 4}, {0, 4}, {0X9999999A, 5}, {0X47AE147B, 5}, {0, 5}, {0X9999999A, 6}, {0X446F8657, 6}, {0, 6},
    {0X970E4F81, 7}, {0X42D6625E, 7}, {0, 7}, {0X970E4F81, 8}, {0X42D6625E, 8}, {0, 8}, {0X966CC01A, 9}, {0X42D6625E, 9}, {0, 9}, {0X966CC01A, 10}, {0X42A38662, 10}, {0, 10},
    {0X966CC01A, 11}, {0X42A38662, 11}, {0, 11}, {0X9E74F884, 12}, {0X4296D1E4, 12}, {0, 12}, {0X9662AB65, 13}, {0X42907805, 13}, {0, 13}, {0X9662AB65, 14}, {0X428D4B2D, 14}, {0, 14},
    {0X9660264C, 15}, {0X428BB4C8, 15}, {0, 15}, {0X9660264C, 16}, {0X428AE996, 16}, {0, 16}, {0X9660264C, 17}, {0X428A83FE, 17}, {0, 17}, {0X9660264C, 18}, {0X428A5132, 18}, {0X2001, 18},
    {0X965FFDFB, 19}, {0X428A37CC, 19}, {0, 19}, {0X965FFDFB, 20}, {0X428A37CC, 20}, {0, 20}, {0X965FF3E7, 21}, {0X428A3172, 21}, {0, 21}, {0X965FEEDD, 22}, {0X428A3172, 22}, {0X201, 22},
    {0X965FEC57, 23}, {0X428A2FDC, 23}, {0, 23}, 
};

static const struct libdivide_u64_branchfree_t denoms64[] = {
    {0, 5}, {0X47AE147AE147AE15, 6}, {0XC71C71C71C71C71D, 7}, {0, 7}, {0X47AE147AE147AE15, 8}, {0XA36E2EB1C432CA58, 9}, {0, 9}, {0X47AE147AE147AE15, 10}, {0XA36E2EB1C432CA58, 11},
    {0, 11}, {0X47AE147AE147AE15, 12}, {0X9B2A7C9FE8B617F1, 13}, {0, 13}, {0X439F40CC28F760CC, 14}, {0X972002FB5C05974D, 15}, {0, 15}, {0X439F40CC28F760CC, 16}, {0X972002FB5C05974D, 17},
    {0, 17}, {0X429E8FC19BD1F6F1, 18}, {0X972002FB5C05974D, 19}, {0, 19}, {0X429E8FC19BD1F6F1, 20}, {0X969FC68151912805, 21}, {0, 21}, {0X429E8FC19BD1F6F1, 22}, {0X969FC68151912805, 23},
    {0, 23}, {0X4F7F449D26A949DD, 24}, {0X967FC0DA51BF3DDD, 25}, {0, 25}, {0X428E8ED5E63B8E8A, 26}, {0X966FBF71EFB24C2F, 27}, {0, 27}, {0X428E8ED5E63B8E8A, 28}, {0X9667BF187DCBD4ED, 29},
    {0, 29}, {0X428A8ECA9A8C6397, 30}, {0X9663BF027395D269, 31}, {0, 31}, {0X428A8ECA9A8C6397, 32}, {0X9661BEFD1A08C80F, 33}, {0, 33}, {0X428A8ECA9A8C6397, 34}, {0X9660BEFBD82195C9, 35},
    {0, 35}, {0X428A8ECA9A8C6397, 36}, {0X96603EFB91E54C07, 37}, {0X40000C000201, 37}, {0X428A4ECA87BD4923, 38}, {0X965FFEFB8574EA53, 39}, {0, 39}, {0X428A4ECA87BD4923, 40},
    {0X965FFEFB8574EA53, 41}, {0, 41}, {0X428A3ECA8603773D, 42}, {0X965FEEFB84B5956C, 43}, {0, 43}, {0X428A36CA8598D950, 44}, {0X965FEEFB84B5956C, 45}, {0X400000C0001, 45},
    {0X428A32CA85801D1A, 46}, {0X965FEAFB84AB8C66, 47}, {0, 47},
};

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
    if (r->cache_bits < 64 && r->stream_size > 0) {
        r->stream_size -= 1;
        uint64_t data = *r->stream++;
#if __LITTLE_ENDIAN__
        data = BSWAP64(data);
#endif
        // Fill in missing bits
        r->data |= data >> r->cache_bits;
        r->next  = data << (64 - r->cache_bits);
        r->cache_bits += 64;
    }
}

static inline uint64_t br_read(br_t* r, size_t num_bits) {
    ASSERT(num_bits <= 64);
    uint64_t shft = 64 - num_bits;
    uint64_t res  = r->data >> shft;

    r->cache_bits -= (uint32_t)num_bits;

    // Append extracted bits from next
    r->data = (r->data << num_bits) | (r->next >> shft);
    r->next <<= num_bits;

    br_load_next(r);

    return res;
}

static inline uint64_t br_peek(br_t* r, size_t num_bits) {
    uint64_t shft = 64 - num_bits;
    uint64_t res  = r->data >> shft;
    return res;
}

static inline br_skip(br_t* r, size_t num_bits) {
    r->cache_bits -= (uint32_t)num_bits;

    // Append extracted bits from next
    r->data = (r->data << num_bits) | (r->next >> (64 - num_bits));
    r->next <<= num_bits;

    br_load_next(r);
}

typedef struct {
    uint64_t size_zy;
    uint32_t size_z;
    uint8_t  num_of_bits;
    uint8_t  part_mask;
    uint8_t  big_shift;
    uint8_t  sml_shift;
    struct libdivide_u64_branchfree_t div_zy;
    struct libdivide_u32_branchfree_t div_z;
} unpack_data_t;

typedef struct {
    union {
        __m128i v;
        int arr[4];
    };
} v4i_t;

static inline v4i_t v4i_set(float x, float y, float z, float w) {
    v4i_t r;
    r.v = _mm_set_epi32(w, z, y, x);
    return r;
}

static inline v4i_t v4i_set1(float x) {
    v4i_t r;
    r.v = _mm_set1_epi32(x);
    return r;
}

static inline v4i_t v4i_add(v4i_t a, v4i_t b) {
    v4i_t r;
    r.v = _mm_add_epi32(a.v, b.v);
    return r;
}

static inline v4i_t v4i_sub(v4i_t a, v4i_t b) {
    v4i_t r;
    r.v = _mm_sub_epi32(a.v, b.v);
    return r;
}

static inline void init_unpack_data_bits(unpack_data_t* data, uint32_t num_of_bits) {
    uint32_t partbits = num_of_bits & 7;
    data->num_of_bits = (uint8_t)num_of_bits;
    data->big_shift = 64 - ((num_of_bits + 7) & ~7);  // Align to the next multiple of 8
    data->sml_shift = (8 - partbits) & 7;
    data->part_mask = partbits ? (1ul << partbits) - 1 : 0xFF;
}

static inline v4i_t unpack_coord(br_t* r, const unpack_data_t* unpack) {
    v4i_t res;
    if (unpack->num_of_bits <= 64) {
        uint64_t w = br_read(r, unpack->num_of_bits);
        uint64_t t = ((w << unpack->sml_shift) & (~0xFFull)) | (w & unpack->part_mask);
        uint64_t v = BSWAP64(t) >> unpack->big_shift;

        uint32_t x1 = libdivide_u64_branchfree_do(v, &unpack->div_zy);
        uint32_t q1 = (uint32_t)(v - x1 * unpack->size_zy);

        uint32_t y1 = libdivide_u32_branchfree_do(q1, &unpack->div_z);
        uint32_t z1 = q1 - y1 * unpack->size_z;

        res.v = _mm_set_epi32(0, z1, y1, x1);
    }

    return res;
}

static inline v4i_t br_read_ints(br_t* br, const uint32_t bitsizes[]) {
    v4i_t res = {0};
    for (int i = 0; i < 3; ++i) {
        _mm_slli_si128(res.v, 32);
        uint32_t num_of_bits = bitsizes[i];
        uint32_t partbits = num_of_bits & 7;
        uint32_t big_shift = 32 - ((num_of_bits + 7) & ~7);  // Align to the next multiple of 8
        uint32_t sml_shift = (8 - partbits) & 7;
        uint32_t part_mask = partbits ? (1ul << partbits) - 1 : 0xFF;
        uint32_t w = (uint32_t)br_read(br, num_of_bits);
        uint32_t r = ((w << sml_shift) & (~0xFFul)) | (w & part_mask);
        uint32_t v = BSWAP32(r) >> big_shift;
        res.v = _mm_blend_epi32(res.v, _mm_cvtsi32_si128(v), 1);
    }
    return res;
} 

static inline write_coord(float* dst, v4i_t coord, float inv_precision) {
    __m128 v = _mm_mul_ps(_mm_cvtepi32_ps(coord.v), _mm_set1_ps(inv_precision));
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);
    _mm_maskstore_ps(dst, mask, v);
}

bool xdr_decompress_coord_float(float* ptr, int* size, float* precision, md_file_o* xfp) {
    int smallidx;
    int natoms;
    int smallnum, smaller, is_smaller, run;
    float *lfp, inv_precision;
    unsigned int bitsize;
    const float* ptrstart = ptr;

    bool result = false;

    v4i_t sizesmall, thiscoord, prevcoord;

    size_t temp_pos = md_temp_get_pos();

    if (xfp == NULL || ptr == NULL) {
        goto done;
    }
    if (!xdr_read_int32(&natoms, 1, xfp)) {
        goto done; /* return if we could not read size */
    }
    *size = natoms;
    uint32_t size3 = natoms * 3;

    /* Dont bother with compression for three atoms or less */
    if (*size <= 9) {
        return xdr_read_float(ptr, size3, xfp);
    }

    /* Compression-time if we got here. Read precision first */
    if (!xdr_read_float(precision, 1, xfp)) {
        goto done;
    }
    v4i_t minint, maxint;
    if (!xdr_read_int32(minint.arr, 3, xfp) || !xdr_read_int32(maxint.arr, 3, xfp)) {
        goto done;
    }

    v4i_t sizeint = v4i_add(v4i_sub(maxint, minint), v4i_set1(1));

    /* check if one of the sizes is to big to be multiplied */
    uint32_t bitsizeint[3] = {0};
    unpack_data_t big_unpack = {0};

    if ((sizeint.arr[0] | sizeint.arr[1] | sizeint.arr[2]) > 0xffffff) {
        bitsizeint[0] = sizeofint(sizeint.arr[0]);
        bitsizeint[1] = sizeofint(sizeint.arr[1]);
        bitsizeint[2] = sizeofint(sizeint.arr[2]);
        bitsize = 0; /* flag the use of large sizes */
    } else {
        bitsize = sizeofints(3, (unsigned int*)sizeint.arr);
        big_unpack.size_zy = sizeint.arr[1] * sizeint.arr[2],
        big_unpack.size_z = sizeint.arr[2],
        big_unpack.num_of_bits = bitsize;
        big_unpack.div_zy = libdivide_u64_branchfree_gen(sizeint.arr[1] * sizeint.arr[2]),
        big_unpack.div_z  = libdivide_u32_branchfree_gen(sizeint.arr[2]),
        init_unpack_data_bits(&big_unpack, bitsize);
    }

    if (!xdr_read_int32(&smallidx, 1, xfp)) {
        return false; /* not sure what has happened here or why we return... */
    }
    int idx = MAX(smallidx - 1, FIRSTIDX);
    smaller = magicints[idx] / 2;
    smallnum = magicints[smallidx] / 2;
    uint32_t smallsize = magicints[smallidx];

    unpack_data_t sml_unpack = {
        .size_zy = (uint64_t)smallsize * smallsize,
        .size_z  = smallsize,
        .div_zy = denoms64[smallidx - FIRSTIDX],
        .div_z  = denoms32[smallidx - FIRSTIDX],
    };

    /* length in bytes */
    uint32_t num_bytes = 0;
    if (!xdr_read_int32((int32_t*)&num_bytes, 1, xfp)) {
        goto done;
    }

    size_t mem_size = ALIGN_TO(num_bytes, 16);
    uint64_t* mem = md_temp_push_aligned(mem_size, 16);

    if (!xdr_read_bytes((uint8_t*)mem, num_bytes, xfp)) {
        goto done;
    }

    br_t br = {0};
    br_init(&br, mem, mem_size / 8);

    lfp = ptr;
    inv_precision = 1.0f / *precision;
    run = 0;
    int atom_idx = 0;
    //printf("bitsize: %i\n", bitsize);
    //printf("sizeint: [%u, %u, %u]\n", sizeint.arr[0], sizeint.arr[1], sizeint.arr[2]);
    int max_smallidx = 0;

    while (atom_idx++ < natoms) {
        if (bitsize == 0) {
            thiscoord = br_read_ints(&br, bitsizeint);
        } else {
            thiscoord = unpack_coord(&br, &big_unpack);
        }

        thiscoord = v4i_add(thiscoord, minint);

        uint64_t data = br_peek(&br, 6);
        uint32_t flag = (uint32_t)(data & 32);
        uint32_t skip = flag ? 6 : 1;
        br_skip(&br, skip);

        is_smaller = 0;
        if (flag) {
            run = (uint32_t)(data & 31);
            is_smaller = run % 3;
            run -= is_smaller;
            is_smaller--;
        }
        if ((lfp - ptrstart) + run > size3) {
            MD_LOG_ERROR("XTC: Buffer overrun during decompression.");
            goto done;
        }
        if (run > 0) {
            //printf("run: %i\n", run);
            //printf("smallidx:  %i\n", smallidx);
            //printf("sizesmall: %i\n", sizesmall.arr[0]);
            //printf("sizesmall: %i\n", sml_unpack.sz);

            prevcoord = thiscoord;

            int count = run / 3;
            v4i_t vsmall = v4i_set1(smallnum);
            v4i_t coord = unpack_coord(&br, &sml_unpack);
            thiscoord = v4i_add(coord, v4i_sub(thiscoord, vsmall));

            write_coord(lfp, thiscoord, inv_precision);
            lfp += 3;
            write_coord(lfp, prevcoord, inv_precision);
            lfp += 3;

            for (int i = 1; i < count; ++i) {
                coord = unpack_coord(&br, &sml_unpack);
                thiscoord = v4i_add(coord, v4i_sub(thiscoord, vsmall));
                write_coord(lfp, thiscoord, inv_precision);
                lfp += 3;
            }
            #if 0
            float ref_x[8], ref_y[8], ref_z[8];

            uint64_t bits[8];
            for (int i = 0; i < count; ++i) {
                bits[i] = br_read_u64(&br, &sml_unpack);
            }

            __m256i v = _mm256_load_si256((__m256i*)bits);
            __m256i x1 = libdivide_u64_branchfree_do_vec256(v, &sml_unpack.div_zy);
            __m256i q1 = _mm256_sub_epi64(v, _mm256_mul_epu32(x1, _mm256_set1_epi32((uint32_t)sml_unpack.size_zy)));
            x1 = _mm256_permutevar8x32_epi32(x1, _mm256_set_epi32(-1, -1, -1, -1, 6, 4, 2, 0));
            q1 = _mm256_permutevar8x32_epi32(q1, _mm256_set_epi32(-1, -1, -1, -1, 6, 4, 2, 0));

            if (count > 4) {
                __m256i w  = _mm256_load_si256((__m256i*)(bits + 4));
                __m256i x2 = libdivide_u64_branchfree_do_vec256(w, &sml_unpack.div_zy);
                __m256i q2 = _mm256_sub_epi64(w, _mm256_mul_epu32(x2, _mm256_set1_epi32((uint32_t)sml_unpack.size_zy)));
                x2 = _mm256_permutevar8x32_epi32(x2, _mm256_set_epi32(6, 4, 2, 0, -1, -1, -1, -1));
                q2 = _mm256_permutevar8x32_epi32(q2, _mm256_set_epi32(6, 4, 2, 0, -1, -1, -1, -1));
                x1 = _mm256_blend_epi32(x1, x2, 0x00AA);
                q1 = _mm256_blend_epi32(q1, q2, 0x00AA);
            }

            __m256i y1 = libdivide_u32_branchfree_do_vec256(q1, &sml_unpack.div_z);
            __m256i z1 = _mm256_sub_epi32(q1, _mm256_mullo_epi32(y1, _mm256_set1_epi32(sml_unpack.size_z)));

            x1 = _mm256_sub_epi32(x1, _mm256_set1_epi32(smallnum));
            y1 = _mm256_sub_epi32(y1, _mm256_set1_epi32(smallnum));
            z1 = _mm256_sub_epi32(z1, _mm256_set1_epi32(smallnum));

            x1 = _mm256_prefix_sum_epi32(x1);
            y1 = _mm256_prefix_sum_epi32(y1);
            z1 = _mm256_prefix_sum_epi32(z1);

            x1 = _mm256_add_epi32(x1, _mm256_set1_epi32(prevcoord.arr[0]));
            y1 = _mm256_add_epi32(y1, _mm256_set1_epi32(prevcoord.arr[1]));
            z1 = _mm256_add_epi32(z1, _mm256_set1_epi32(prevcoord.arr[2]));

            __m256 scl = _mm256_set1_ps(inv_precision);
            __m256 x = _mm256_mul_ps(_mm256_cvtepi32_ps(x1), scl);
            __m256 y = _mm256_mul_ps(_mm256_cvtepi32_ps(y1), scl);
            __m256 z = _mm256_mul_ps(_mm256_cvtepi32_ps(z1), scl);

            // Interleave x, y, z and write back to memory

            // Step 1: Interleave x and y using unpack operations
            __m256 xy_lo = _mm256_unpacklo_ps(x, y); // [ x0, y0, x1, y1, x2, y2, x3, y3 ]
            __m256 xy_hi = _mm256_unpackhi_ps(x, y); // [ x4, y4, x5, y5, x6, y6, x7, y7 ]

            // Step 2: Interleave z into the xy pairs
            __m256 xyz_lo = _mm256_unpacklo_ps(xy_lo, z); // [ x0, y0, z0, x1, y1, z1, x2, y2 ]
            __m256 xyz_hi = _mm256_unpackhi_ps(xy_lo, z); // [ z2, x3, y3, z3, x4, y4, z4, x5 ]
            __m256 xyz_ex = _mm256_unpackhi_ps(xy_hi, z); // [ y5, z5, x6, y6, z6, x7, y7, z7 ]

            // Store first entry
            _mm_storeu_ps(lfp, _mm256_castps256_ps128(xyz_lo));
            lfp += 3;

            // Swap in prevcoords XYZ
            __m128 prev_xyz = _mm_mul_ps(_mm_cvtepi32_ps(prevcoord.v), _mm256_castps256_ps128(scl));
            xyz_lo = _mm256_blend_ps(xyz_lo, _mm256_castps128_ps256(prev_xyz), 7);

            // Store count elements to
            _mm256_storeu_ps(lfp, xyz_lo);
            _mm256_storeu_ps(lfp + 8, xyz_hi);
            _mm256_storeu_ps(lfp + 16, xyz_ex);
            lfp += run;
            #endif
            atom_idx += count;
        } else {
            write_coord(lfp, thiscoord, inv_precision);
            lfp += 3;
        }
        smallidx += is_smaller;
        if (is_smaller < 0) {
            smallnum = smaller;
            smaller = (smallidx > FIRSTIDX) ? magicints[smallidx - 1] / 2 : 0;
        } else if (is_smaller > 0) {
            smaller = smallnum;
            smallnum = magicints[smallidx] / 2;
        }
        uint32_t new_small = magicints[smallidx];
        if (new_small == 0) {
            MD_LOG_ERROR("XTC: Invalid size found in 'xdrfile_decompress_coord_float'.");
            goto done;
        }
        if (new_small != sml_unpack.size_z) {
            uint32_t idx = smallidx - FIRSTIDX;
            sml_unpack.size_zy      = (uint64_t)new_small * new_small;
            sml_unpack.size_z       = new_small;
            sml_unpack.num_of_bits  = smallidx;
            sml_unpack.div_zy       = denoms64[idx];
            sml_unpack.div_z        = denoms32[idx];

            init_unpack_data_bits(&sml_unpack, smallidx);
        }
    }
    result = true;
done:
    md_temp_set_pos_back(temp_pos);
    return result;
}




static bool xtc_frame_header(XDRFILE* xd, int* natoms, int* step, float* time, float box[3][3]) {
    int magic;

    if ((xdrfile_read_int(&magic, 1, xd)) != 1) {
        MD_LOG_ERROR("XTC: Failed to read magic number in header");
        return false;
    }
    if (magic != XTC_MAGIC) {
        MD_LOG_ERROR("XTC: Magic number did not match");
        return false;
    }
    if ((xdrfile_read_int(natoms, 1, xd)) != 1) {
        MD_LOG_ERROR("XTC: Failed to read number of atoms");
        return false;
    }
    if ((xdrfile_read_int(step, 1, xd)) != 1) {
        MD_LOG_ERROR("XTC: Failed to read step");
        return false;
    }
    if ((xdrfile_read_float(time, 1, xd)) != 1) {
        MD_LOG_ERROR("XTC: Failed to read timestamp");
        return false;
    }
    if ((xdrfile_read_float((float*)box, DIM * DIM, xd)) != DIM * DIM) {
        MD_LOG_ERROR("XTC: Failed to read box dimensions");
        return false;
    }

    return true;
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

bool md_xtc_read_frame_offsets_and_times(md_file_o* xdr, md_array(int64_t)* frame_offsets, md_array(double)* frame_times, md_allocator_i* alloc) {
    int step, natoms;
    float time;
    float box[3][3];

    size_t filesize = (size_t)md_file_size(xdr);

    if (filesize == 0) {
        MD_LOG_ERROR("XTC: Failed extract filesize");
        return false;
    }

    /* Go to file beg */
    if (!md_file_seek(xdr, 0, MD_FILE_BEG)) {
        MD_LOG_ERROR("XTC: Failed to seek to beginning of file");
        return false;
    }

    if (!md_xtc_read_frame_header(xdr, &natoms, &step, &time, box)) {
        return false;
    }

    /* Dont bother with compression for nine atoms or less */
    if (natoms <= 9) {
        const size_t framebytes = XTC_SMALL_HEADER_SIZE + XTC_SMALL_COORDS_SIZE * natoms;
        const size_t num_frames = (filesize / framebytes); /* Should we complain if framesize doesn't divide filesize? */

        md_array_ensure(*frame_offsets, num_frames, alloc);
        md_array_ensure(*frame_times,   num_frames, alloc);

        md_array_push(*frame_offsets, 0, alloc);
        md_array_push(*frame_times, time, alloc);

        for (size_t i = 1; i < num_frames; i++) {
            const size_t offset = i * framebytes;

            if (!md_file_seek(xdr, offset, MD_FILE_BEG) || !md_xtc_read_frame_header(xdr, &natoms, &step, &time, box)) {
                goto fail;
            }

            md_array_push(*frame_offsets, i * framebytes, alloc);
            md_array_push(*frame_times, time, alloc);
        }
        md_array_push(*frame_offsets, num_frames * framebytes, alloc);

        return true;
    } else {
        int framebytes, est_nframes;

        /* Move pos back to end of first header */
        if (!md_file_seek(xdr, (int64_t)XTC_HEADER_SIZE, MD_FILE_BEG)) {
            return false;
        }

        if (!xdr_read_int32(&framebytes, 1, xdr)) {
            MD_LOG_ERROR("XTC: Failed to read framebytes");
            return false;
        }
        framebytes = (framebytes + 3) & ~0x03; /* Rounding to the next 32-bit boundary */
        est_nframes = (int)(filesize / ((int64_t)(framebytes + XTC_HEADER_SIZE)) + 1); /* must be at least 1 for successful growth */
        /* First `framebytes` might be larger than average, so we would underestimate `est_nframes`
        */
        est_nframes += est_nframes / 5;

        /* Skip `framebytes` */
        if (!md_file_seek(xdr, (int64_t)(framebytes), MD_FILE_CUR)) {
            goto fail;
        }

        md_array_ensure(*frame_offsets, (size_t)est_nframes, alloc);
        md_array_ensure(*frame_times,   (size_t)est_nframes, alloc);

        md_array_push(*frame_offsets, 0, alloc);
        md_array_push(*frame_times, time, alloc);

        while (true) {
            const int64_t offset = md_file_tell(xdr);
            if (offset == (int64_t)filesize) {
                // Add last offset
                md_array_push(*frame_offsets, offset, alloc);
                return true;
            }
            if (!md_xtc_read_frame_header(xdr, &natoms, &step, &time, box)) {
                goto fail;
            }

            /* Store position in `offsets`, adjust for header */
            md_array_push(*frame_offsets, offset, alloc);
            md_array_push(*frame_times, time, alloc);

            if (!md_file_seek(xdr, offset + XTC_HEADER_SIZE, MD_FILE_BEG)) {
                goto fail;
            }
            /* Read how much to skip */
            if (!xdr_read_int32(&framebytes, 1, xdr)) {
                goto fail;
            }
            
            framebytes = (framebytes + 3) & ~0x03; /* Rounding to the next 32-bit boundary */

            /* Skip `framebytes` to next header */
            if (!md_file_seek(xdr, framebytes, MD_FILE_CUR)) {
                goto fail;
            }
        }

    fail:
        md_array_free(*frame_offsets, alloc);
        md_array_free(*frame_times, alloc);
        *frame_offsets = 0;
        *frame_times = 0;
        return false;
    }
}

static bool xtc_frame_coords(XDRFILE* xd, int natoms, rvec* x) {
    int result;
    float prec;

    result = xdrfile_decompress_coord_float(x[0], &natoms, &prec, xd);
    if (result != natoms) {
        MD_LOG_ERROR("XTC: Failed to read coordinates");
        return false;
    }

    return true;
}

size_t md_xtc_read_frame_coords(md_file_o* xdr_file, float* coords, size_t capacity) {
    int natoms;
    float prec;

    bool result = xdr_decompress_coord_float(coords, &natoms, &prec, xdr_file);
    if (!result) {
        MD_LOG_ERROR("XTC: Failed to read coordinates");
        return false;
    }

    return true;
}

static size_t xtc_frame_offsets_and_times(XDRFILE* xd, md_array(int64_t)* frame_offsets, md_array(double)* frame_times, md_allocator_i* alloc) {
    int step, natoms;
    float time;
    float box[3][3];

    /* Go to file beg */
    if (xdr_seek(xd, 0L, SEEK_SET) != exdrOK) {
        MD_LOG_ERROR("XTC: Failed to seek to beginning of file");
        return 0;
    }

    if (!xtc_frame_header(xd, &natoms, &step, &time, box)) {
        return 0;
    }

    /* Go to file end */
    if (xdr_seek(xd, 0L, SEEK_END) != exdrOK) {
        return 0;
    }
    /* Cursor position is equivalent to file size */
    size_t filesize = (size_t)xdr_tell(xd);

    /* Dont bother with compression for nine atoms or less */
    if (natoms <= 9) {
        const size_t framebytes = XTC_SMALL_HEADER_SIZE + XTC_SMALL_COORDS_SIZE * natoms;
        const size_t num_frames = (filesize / framebytes); /* Should we complain if framesize doesn't divide filesize? */

        md_array_push(*frame_offsets, 0, alloc);
        md_array_push(*frame_times, time, alloc);

        for (size_t i = 1; i < num_frames; i++) {
            const size_t offset = i * framebytes;

            if (xdr_seek(xd, offset, SEEK_SET) != exdrOK || !xtc_frame_header(xd, &natoms, &step, &time, box)) {
                md_array_free(*frame_offsets, alloc);
                md_array_free(*frame_times, alloc);
                *frame_offsets = 0;
                *frame_times = 0;
                return 0;
            }

            md_array_push(*frame_offsets, i * framebytes, alloc);
            md_array_push(*frame_times, time, alloc);
        }
        md_array_push(*frame_offsets, num_frames * framebytes, alloc);

        return num_frames;
    } else {
        int framebytes, est_nframes;

        /* Move pos back to end of first header */
        if (xdr_seek(xd, (int64_t)XTC_HEADER_SIZE, SEEK_SET) != exdrOK) {
            return false;
        }
        
        if (xdrfile_read_int(&framebytes, 1, xd) == 0) {
            MD_LOG_ERROR("XTC: Failed to read framebytes");
            return false;
        }
        framebytes = (framebytes + 3) & ~0x03; /* Rounding to the next 32-bit boundary */
        est_nframes = (int)(filesize / ((int64_t)(framebytes + XTC_HEADER_SIZE)) +
            1); /* must be at least 1 for successful growth */
                /* First `framebytes` might be larger than average, so we would underestimate `est_nframes`
                */
        est_nframes += est_nframes / 5;

        /* Skip `framebytes` */
        if (xdr_seek(xd, (int64_t)(framebytes), SEEK_CUR) != exdrOK) {
            goto fail;
        }
        
        md_array_ensure(*frame_offsets, (size_t)est_nframes, alloc);
        md_array_ensure(*frame_times,   (size_t)est_nframes, alloc);

        md_array_push(*frame_offsets, 0, alloc);
        md_array_push(*frame_times, time, alloc);
        size_t num_frames = 1;

        while (1) {
            const int64_t offset = xdr_tell(xd);
            if (offset == (int64_t)filesize) {
                // Add last offset
                md_array_push(*frame_offsets, offset, alloc);
                return num_frames;
            }
            if (!xtc_frame_header(xd, &natoms, &step, &time, box)) {
                goto fail;
            }

            /* Store position in `offsets`, adjust for header */
            md_array_push(*frame_offsets, offset, alloc);
            md_array_push(*frame_times, time, alloc);
            num_frames += 1;

            if (xdr_seek(xd, offset + XTC_HEADER_SIZE, SEEK_SET) != exdrOK) {
                goto fail;
            }
            /* Read how much to skip */
            if (xdrfile_read_int(&framebytes, 1, xd) == 0) {
                goto fail;
            }
            framebytes = (framebytes + 3) & ~0x03; /* Rounding to the next 32-bit boundary */

            /* Skip `framebytes` to next header */
            if (xdr_seek(xd, framebytes, SEEK_CUR) != exdrOK) {
                goto fail;
            }
        }

    fail:
        md_array_free(*frame_offsets, alloc);
        md_array_free(*frame_times, alloc);
        *frame_offsets = 0;
        *frame_times = 0;
        return 0;
    }
}

bool xtc_get_header(struct md_trajectory_o* inst, md_trajectory_header_t* header) {
    xtc_t* xtc = (xtc_t*)inst;
    ASSERT(xtc);
    ASSERT(xtc->magic == MD_XTC_TRAJ_MAGIC);
    ASSERT(header);

    *header = xtc->header;
    return true;
}

// This is lowlevel cruft for enabling parallel loading and decoding of frames
static size_t xtc_fetch_frame_data(struct md_trajectory_o* inst, int64_t frame_idx, void* frame_data_ptr) {
    xtc_t* xtc = (xtc_t*)inst;
    ASSERT(xtc);
    ASSERT(xtc->magic == MD_XTC_TRAJ_MAGIC);

    if (!xtc->file) {
        MD_LOG_ERROR("XTC: File handle is NULL");
        return 0;
    }

    if (!xtc->frame_offsets) {
        MD_LOG_ERROR("XTC: Frame offsets is empty");
        return 0;
    }

    if (frame_idx < 0 || (int64_t)xtc->header.num_frames <= frame_idx) {
        MD_LOG_ERROR("XTC: Frame index is out of range");
        return 0;
    }

    const int64_t beg = xtc->frame_offsets[frame_idx];
    const int64_t end = xtc->frame_offsets[frame_idx + 1];
    const size_t frame_size = (size_t)MAX(0, end - beg);

    if (frame_data_ptr) {
        ASSERT(xtc->file);
        md_mutex_lock(&xtc->mutex);
        // Seek and read must be an atomic operation to avoid race conditions
        xdr_seek(xtc->file, beg, SEEK_SET);
        const size_t bytes_read = xdr_read(xtc->file, frame_data_ptr, frame_size);
        md_mutex_unlock(&xtc->mutex);
        (void)bytes_read;
        ASSERT(frame_size == bytes_read);
    }
    return frame_size;
}

static bool xtc_decode_frame_data(struct md_trajectory_o* inst, const void* frame_data_ptr, size_t frame_data_size, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(inst);
    ASSERT(frame_data_ptr);
    ASSERT(frame_data_size);

    bool result = true;

    xtc_t* xtc = (xtc_t*)inst;
    if (xtc->magic != MD_XTC_TRAJ_MAGIC) {
        MD_LOG_ERROR("XTC: Error when decoding frame coord, xtc magic did not match");
        return false;
    }

    if ((x || y || z) && !(x && y && z)) {
        MD_LOG_ERROR("XTC: User supplied coordinates (x,y,z) cannot be partially supplied");
        return false;
    }

    // There is a warning for ignoring const qualifier for frame_data_ptr, but it is ok since we only perform read operations "r" with the data.
    XDRFILE* file = xdrfile_mem((void*)frame_data_ptr, frame_data_size, "r");
    ASSERT(file);

    // Get header
    int natoms = 0, step = 0;
    float time = 0;
    float box[3][3];
    result = xtc_frame_header(file, &natoms, &step, &time, box);
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
            size_t byte_size = natoms * sizeof(rvec);
            rvec* pos = md_alloc(md_get_heap_allocator(), byte_size);
            result = xtc_frame_coords(file, natoms, pos);
            if (result) {            
                // nm -> Ångström
                for (int i = 0; i < natoms; ++i) {
                    x[i] = pos[i][0] * 10.0f;
                    y[i] = pos[i][1] * 10.0f;
                    z[i] = pos[i][2] * 10.0f;
                }
            }
            md_free(md_get_heap_allocator(), pos, byte_size);
        }
    }

    xdrfile_close(file);

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
    const size_t frame_size = xtc_fetch_frame_data(inst, frame_idx, NULL);
    if (frame_size > 0) {
        md_allocator_i* alloc = md_get_heap_allocator();

        void* frame_data = md_alloc(alloc, frame_size);
        ASSERT(frame_data);

        const size_t read_size = xtc_fetch_frame_data(inst, frame_idx, frame_data);
        if (read_size != frame_size) {
            MD_LOG_ERROR("Failed to read the expected size");
            goto done;
        }
        result = xtc_decode_frame_data(inst, frame_data, frame_size, header, x, y, z);
    done:
        md_free(alloc, frame_data, frame_size);
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

    // Ensure that the path is zero terminated (not guaranteed by str_t)
    md_strb_t sb = md_strb_create(md_get_temp_allocator());
    md_strb_push_str(&sb, filename);
    XDRFILE* file = xdrfile_open(md_strb_to_cstr(sb), "r");

    if (file) {
        xdr_seek(file, 0L, SEEK_END);
        const size_t filesize = (size_t)xdr_tell(file);
        xdr_seek(file, 0L, SEEK_SET);

        int num_atoms, step;
        float time;
        float box[3][3];
        if (!xtc_frame_header(file, &num_atoms, &step, &time, box)) {
            goto fail;
        }

        if (num_atoms == 0) {
            MD_LOG_ERROR("XTC: Number of atoms in trajectory was zero");
            goto fail;
        }

        md_strb_push_str(&sb, STR_LIT(".cache"));
        str_t path = md_strb_to_str(sb);

        xtc_cache_t cache = {0};
        if (!try_read_cache(&cache, path, filesize, alloc)) {
            cache.header.magic     = MD_XTC_CACHE_MAGIC;
            cache.header.version   = MD_XTC_CACHE_VERSION;
            cache.header.num_bytes = filesize;
            cache.header.num_atoms = num_atoms;
            cache.header.num_frames = xtc_frame_offsets_and_times(file, &cache.frame_offsets, &cache.frame_times, alloc);
            if (!cache.header.num_frames) {
                goto fail;
            }
            if (!cache.frame_offsets || !cache.frame_times) {
                MD_LOG_DEBUG("XTC: frame offsets or frame times was empty");
                goto fail;
            }

            if (!(flags & MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE)) {
                // If we fail to write the cache, that's ok, we can inform about it, but do not halt
                if (write_cache(&cache, path)) {
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
        xtc->frame_offsets = cache.frame_offsets;
        xtc->mutex = md_mutex_create();

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
        //traj->fetch_frame_data = xtc_fetch_frame_data;
        //traj->decode_frame_data = xtc_decode_frame_data;

        return traj;
    }
fail:
    if (file) xdrfile_close(file);
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

    if (xtc->file) xdrfile_close(xtc->file);
    md_mutex_destroy(&xtc->mutex);
    md_arena_allocator_destroy(xtc->allocator);
}

static md_trajectory_loader_i xtc_loader = {
    md_xtc_trajectory_create,
    md_xtc_trajectory_free,
};

md_trajectory_loader_i* md_xtc_trajectory_loader(void) {
    return &xtc_loader;
}

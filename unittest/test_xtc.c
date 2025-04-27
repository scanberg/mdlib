#include "utest.h"
#include <string.h>

#include <md_xtc.h>
#include <md_trajectory.h>
#include <md_molecule.h>
#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_log.h>

#include <xdrfile_xtc.h>

#define FULL_TEST 0

static inline uint64_t decodebits(int buf[3], int num_of_bits) {
    int cnt;
    unsigned char* cbuf;
    uint64_t lastbits, lastbyte;
    uint64_t mask = num_of_bits < 64 ? (1LLU << num_of_bits) - 1 : 0xFFFFFFFFFFFFFFFFLLU;

    cbuf = ((unsigned char*)buf) + 3 * sizeof(*buf);
    cnt = buf[0];
    lastbits = buf[1];
    lastbyte = buf[2];

    uint64_t num = 0;
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

typedef struct br_t {
    const uint64_t* stream;
    uint64_t data;
    uint64_t next;
    uint32_t cache_bits;
    uint32_t stream_size;
} br_t;

static void br_init(br_t* r, const uint64_t* stream, size_t num_qwords) {
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

static inline uint64_t br_peek(br_t* r, size_t num_bits) {
    ASSERT(num_bits <= 64);
    uint64_t shft = 64 - num_bits;
    uint64_t res  = r->data >> shft;
    return res;
}

static inline void br_load_next(br_t* r) {
#if 0
    if (r->cache_bits <= 64 && r->stream_size > 0) {
        r->stream_size -= 1;
        uint64_t data = *r->stream++;
#if __LITTLE_ENDIAN__
        data = BSWAP64(data);
#endif
        if (r->cache_bits < 64) {
            // Fill in missing bits
            r->data |= data >> r->cache_bits;
            r->next  = data << (64 - r->cache_bits);
        } else {
            r->next = data;
        }
        r->cache_bits += 64;
    }
#else
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
#endif
}

static inline uint64_t br_read(br_t* r, size_t num_bits) {
#if 0
    ASSERT(num_bits <= 64);
    uint64_t shft = 64 - num_bits;
    uint64_t res  = r->data >> shft;

    r->cache_bits -= (uint32_t)num_bits;

    // Append extracted bits from next
    if (num_bits < 64) {
        r->data = (r->data << num_bits) | (r->next >> shft);
        r->next <<= num_bits;
    } else {
        r->data = r->next;
        r->next = 0;
    }

    if (r->cache_bits <= 64 && r->stream_size > 0) {
        r->stream_size -= 1;
        uint64_t data = *r->stream++;
#if __LITTLE_ENDIAN__
        data = BSWAP64(data);
#endif
        // Fill in missing bits
        if (r->cache_bits < 64) {
            r->data |= data >> r->cache_bits;
            r->next  = data << (64 - r->cache_bits);
        } else {;
            r->next  = data;
        }
        r->cache_bits += 64;
    }

    return res;
#else
    ASSERT(num_bits <= 64);

    uint64_t res = r->data >> (64 - num_bits);
    r->cache_bits -= (uint32_t)num_bits;

    // Unified handling
    if (num_bits < 64) {
        r->data <<= num_bits;
        r->data |= r->next >> (64 - num_bits);
        r->next <<= num_bits;
    } else {
        r->data = r->next;
        //r->next = 0;
    }

    br_load_next(r);
    return res;
#endif
}

static const int num_bits[] = {5, 56, 48, 7, 1, 2, 64, 5, 32, 8, 55, 8, 55, 55, 48, 8, 4, 1, 1, 2, 3, 7, 64, 64, 64, 64, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32 ,33, 35, 41, 44, 51, 55, 59, 63, 63, 63, 1, 2, 3,4};

static const int num_bits_32[] = {5, 26, 28, 7, 1, 2, 14, 5, 32, 8, 15, 8, 15, 25, 18, 8, 4, 1, 1, 2, 3, 7, 32, 32, 32, 32, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32 ,23, 25, 21, 24, 21, 25, 29, 23, 23, 13, 1, 2, 3,4};


UTEST(xtc, bitread) {
    size_t temp_pos = md_temp_get_pos();
    int* buf = md_temp_push(sizeof(int) * (1024 + 3));
    buf[0] = buf[1] = buf[2] = 0;
    srand(0);
    for (int i = 3; i < 1024; ++i) {
        buf[i] = rand() << 16 | rand();
    }

    br_t r;
    br_init(&r, (const uint64_t*)&buf[3], 512 * sizeof(int));

    size_t bit_offset = 0;
    uint8_t* base = &buf[3];

    for (size_t i = 0; i < ARRAY_SIZE(num_bits_32); ++i) {
        size_t num_of_bits = num_bits_32[i];
        uint64_t ref = decodebits(buf, num_of_bits);
        uint64_t val = br_read(&r, num_of_bits);
        EXPECT_EQ(val, ref);

        uint64_t byte_offset = bit_offset >> 3;
        uint64_t shift       = bit_offset & 7;

        uint64_t x;
        MEMCPY(&x, base + byte_offset, 8);
        x = BSWAP64(x);
        x <<= shift;
        x >>= (64 - num_of_bits);

        bit_offset += num_of_bits;
    }

    for (uint32_t num_of_bits = 1; num_of_bits <= 64; ++num_of_bits) {
        buf[0] = buf[1] = buf[2] = 0;
        br_init(&r, (const uint64_t*)&buf[3], 512 * sizeof(int));
        for (size_t i = 0; i < 512; ++i) {
            uint64_t ref = decodebits(buf, num_of_bits);
            uint64_t val = br_read(&r,     num_of_bits);
            EXPECT_EQ(val, ref);
            if (val != ref) {
                printf("Error with N=%i\n", num_of_bits);
            }
        }
    }

    md_temp_set_pos_back(temp_pos);
}

__m256i bswap32_avx2(__m256i x) {
    const __m256i shuffle_mask = _mm256_set_epi8(
        12,13,14,15,  8, 9,10,11,
        4, 5, 6, 7,  0, 1, 2, 3,
        12,13,14,15,  8, 9,10,11,
        4, 5, 6, 7,  0, 1, 2, 3
    );
    return _mm256_shuffle_epi8(x, shuffle_mask);
}

__m256i bswap64_avx2(__m256i x) {
    const __m256i shuffle_mask = _mm256_set_epi8(
        8, 9,10,11,12,13,14,15,
        0, 1, 2, 3, 4, 5, 6, 7,
        8, 9,10,11,12,13,14,15,
        0, 1, 2, 3, 4, 5, 6, 7
    );
    return _mm256_shuffle_epi8(x, shuffle_mask);
}

// Only accurate up to 56 bits
__m256i load_packed_bits_avx2_4_56_be(const uint8_t* data, int bit_offset, int N) {
    __m256i big_shift = md_mm256_set1_epi64(64 - N);

    // Compute bit offsets per lane
    __m256i idx         = md_mm256_set_epi64(3,2,1,0);
    __m256i N_vec       = md_mm256_set1_epi64(N);
    __m256i base_bit    = md_mm256_set1_epi64(bit_offset);
    __m256i bit_offsets = _mm256_add_epi64(base_bit, _mm256_mul_epu32(idx, N_vec));

    // Byte offset for each lane
    __m256i byte_offsets = _mm256_srli_epi64(bit_offsets, 3);

    // Bit shifts within each 32-bit load
    __m256i shift      = _mm256_and_si256(bit_offsets, md_mm256_set1_epi64(7));

    // Load 8 bytes from base
    __m256i gather     = _mm256_i64gather_epi64((const __int64_t*)data, byte_offsets, 1);

    __m256i shifted    = _mm256_sllv_epi64(gather, shift);
    return _mm256_srlv_epi64(shifted, big_shift);
}

typedef struct {
    uint64_t lo_mask;
    uint64_t hi_mask;
    int lo_shift;
    int hi_shift;
    int big_shift;
    int num_of_bits;
} unpack_t;

__m256i load_packed_bits_avx2_4_56_le(const uint8_t* data, int bit_offset, const unpack_t* unpack) {
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

    __m256i y = bswap64_avx2(x);

    __m256i hi_shift = _mm256_sub_epi64(_mm256_set1_epi64x(unpack->hi_shift), shift);
    __m256i lo_shift = _mm256_sub_epi64(_mm256_set1_epi64x(unpack->lo_shift), shift);

    __m256i hi_mask  = _mm256_set1_epi64x(unpack->hi_mask);
    __m256i lo_mask  = _mm256_set1_epi64x(unpack->lo_mask);

    __m256i hi = _mm256_and_si256(_mm256_srlv_epi64(y, hi_shift), hi_mask);
    __m256i lo = _mm256_and_si256(_mm256_srlv_epi64(y, lo_shift), lo_mask);

    __m256i combined = _mm256_or_si256(hi, lo);

    __m256i swapped  = bswap64_avx2(combined);

    __m256i result   = _mm256_srlv_epi64(swapped, md_mm256_set1_epi64(unpack->big_shift));
    return result;
}

UTEST(xtc, decode_bits) {
    size_t temp_pos = md_temp_get_pos();
    int* buf = md_temp_push(sizeof(int) * (1024 + 3));
    buf[0] = buf[1] = buf[2] = 0;
    srand(0);
    for (int i = 3; i < 1024; ++i) {
        buf[i] = rand() << 16 | rand();
    }

    br_t r = {0};
    br_init(&r, (const uint64_t*)&buf[3], 512 * sizeof(int));
    
    size_t bit_offset = 0;
    uint8_t* base = &buf[3];

    for (size_t i = 0; i < ARRAY_SIZE(num_bits_32); ++i) {
        int num_of_bits = num_bits_32[i];
           
        int fullbytes = num_of_bits >> 3;
        int partbits  = num_of_bits &  7;

        uint64_t big_shift = 64 - ((num_of_bits + 7) & ~7);  // Align to the next multiple of 8
        uint64_t sml_shift = (8 - partbits) & 7;             // Avoid branch, ensures zero when partbits == 0
        uint64_t part_mask = partbits ? (1ull << partbits) - 1 : 0xFF;

        unpack_t unpack = {
            .lo_mask     = part_mask,
            .hi_mask     = ((1LLU << (num_of_bits + sml_shift)) - 1) & (~0xFFLLU),
            .lo_shift    = 64 - num_of_bits,
            .hi_shift    = 64 - num_of_bits - sml_shift,
            .big_shift   = big_shift,
            .num_of_bits = num_of_bits,
        };


        uint64_t v = 0;
        int i = 0;
        for (; i < fullbytes; i++) {
            uint64_t ibyte = decodebits(buf, 8);
            v |= ibyte << (8 * i);
        }

        if (partbits) {
            v |= ((uint64_t) decodebits(buf, partbits)) << (8 * i);
        }



        uint64_t byte_offset = bit_offset >> 3;
        uint64_t shift       = bit_offset  & 7;

        uint64_t raw;
        MEMCPY(&raw, base + byte_offset, 8);

        uint64_t x = BSWAP64(raw);

        uint32_t z = BSWAP32(raw);
        z <<= shift;
        z >>= (32 - num_of_bits);

        uint64_t y = x;


        x <<= shift;
        x >>= (64 - num_of_bits);

        int hi_shift = 64 - num_of_bits - shift - sml_shift;
        int lo_shift = 64 - num_of_bits - shift;

        uint64_t hi_mask = ((1LLU << (num_of_bits + sml_shift)) - 1) & (~0xFFLLU);
        uint64_t lo_mask = part_mask;

        uint64_t u = (y >> hi_shift) & hi_mask | ((y >> lo_shift) & lo_mask);
        uint64_t e = BSWAP64(u) >> (big_shift);

        // Remove unnecessary masking by directly applying shifts
        uint64_t t = (x << sml_shift) & ~0xFFull | (x & part_mask);
        // Byte-swap and shift to align the result correctly
        uint64_t q = BSWAP64(t) >> big_shift;

        EXPECT_EQ(v, q);
        if (v != q) {
            printf("Error with N=%i\n", num_of_bits);
        }

        __m256i cool = load_packed_bits_avx2_4_56_le(base, bit_offset, &unpack);

        bit_offset += num_of_bits;
    }

    md_temp_set_pos_back(temp_pos);
}

UTEST(xtc, trajectory_i) {
    const str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.xtc");
    md_trajectory_i* traj = md_xtc_trajectory_create(path, md_get_heap_allocator(), MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE);
    ASSERT_TRUE(traj);

    const int64_t num_atoms  = md_trajectory_num_atoms(traj);
    const int64_t num_frames = md_trajectory_num_frames(traj);

    EXPECT_EQ(num_atoms, 1336);
    EXPECT_EQ(num_frames, 501);

    const int64_t mem_size = num_atoms * 3 * sizeof(float);
    void* mem_ptr = md_alloc(md_get_heap_allocator(), mem_size);
    float *x = (float*)mem_ptr;
    float *y = (float*)mem_ptr + num_atoms * 1;
    float *z = (float*)mem_ptr + num_atoms * 2;

    md_trajectory_frame_header_t header;

    for (int64_t i = 0; i < num_frames; ++i) {
        EXPECT_TRUE(md_trajectory_load_frame(traj, i, &header, x, y, z));
    }

    md_free(md_get_heap_allocator(), mem_ptr, mem_size);
    md_xtc_trajectory_free(traj);
}

UTEST(xtc, catalyst) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));

    const str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.xtc");
    md_file_o*  file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);

    XDRFILE* xdr = xdrfile_open(path.ptr, "r");

    ASSERT_TRUE(file);
    ASSERT_TRUE(xdr);

    md_array(int64_t) frame_offsets = 0;
    md_array(double)  frame_times = 0;

    const size_t num_atoms  = 1336;
    const size_t num_frames = 501;

    const size_t coord_size = num_atoms * 3 * sizeof(float);
    float *ref = (float*)md_vm_arena_push(arena, coord_size);
    float *xyz = (float*)md_vm_arena_push(arena, coord_size);

    md_xtc_read_frame_offsets_and_times(file, &frame_offsets, &frame_times, arena);
    size_t xtc_num_frames = frame_offsets ? md_array_size(frame_offsets) - 1 : 0;
    EXPECT_EQ(xtc_num_frames, num_frames);

    int natoms, step;
    float time, box[3][3];

    for (size_t i = 0; i < num_frames; ++i) {
        md_file_seek(file, frame_offsets[i], MD_FILE_BEG);
        EXPECT_TRUE(md_xtc_read_frame_header(file, &natoms, &step, &time, box));
        EXPECT_TRUE(md_xtc_read_frame_coords(file, (float*)xyz, num_atoms));

        static const size_t xtc_header_size = 52;
        xdr_seek(xdr, frame_offsets[i] + xtc_header_size, SEEK_SET);
        int ncoord = natoms;
        float prec;
        if (xdrfile_decompress_coord_float(ref, &ncoord, &prec, xdr) != natoms) {
            MD_LOG_ERROR("Error reading coordinates from XDR file\n");
            goto done;
        }
        EXPECT_EQ(ncoord, num_atoms);

        for (int j = 0; j < num_atoms; ++j) {
            EXPECT_EQ(ref[j * 3 + 0], xyz[j * 3 + 0]);
            EXPECT_EQ(ref[j * 3 + 1], xyz[j * 3 + 1]);
            EXPECT_EQ(ref[j * 3 + 2], xyz[j * 3 + 2]);
        }
    }
    
done:
    md_vm_arena_destroy(arena);
    md_file_close(file);
    xdrfile_close(xdr);
}

UTEST(xtc, amyloid) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));

    const str_t path = STR_LIT("E:/data/md/amyloid-6T/prod-centered.xtc");
    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);

    XDRFILE* xdr = xdrfile_open(path.ptr, "r");

    ASSERT_TRUE(file);
    ASSERT_TRUE(xdr);

    md_array(int64_t) frame_offsets = 0;
    md_array(double)  frame_times = 0;

    const size_t num_atoms  = 161271;
    const size_t num_frames = 5701;

    const size_t coord_size = num_atoms * 3 * sizeof(float);
    float *ref = (float*)md_vm_arena_push(arena, coord_size);
    float *xyz = (float*)md_vm_arena_push(arena, coord_size);

    md_xtc_read_frame_offsets_and_times(file, &frame_offsets, &frame_times, arena);
    size_t xtc_num_frames = frame_offsets ? md_array_size(frame_offsets) - 1 : 0;
    EXPECT_EQ(xtc_num_frames, num_frames);

    int natoms, step;
    float time, box[3][3];

    for (size_t i = 0; i < 5; ++i) {
        md_file_seek(file, frame_offsets[i], MD_FILE_BEG);
        EXPECT_TRUE(md_xtc_read_frame_header(file, &natoms, &step, &time, box));
        EXPECT_TRUE(md_xtc_read_frame_coords(file, xyz, num_atoms));

        static const size_t xtc_header_size = 52;
        xdr_seek(xdr, frame_offsets[i] + xtc_header_size, SEEK_SET);
        int ncoord = natoms;
        float prec;
        if (xdrfile_decompress_coord_float(ref, &ncoord, &prec, xdr) != natoms) {
            MD_LOG_ERROR("Error reading coordinates from XDR file\n");
            goto done;
        }
        EXPECT_EQ(ncoord, num_atoms);

        for (int j = 0; j < num_atoms; ++j) {
            EXPECT_EQ(ref[j * 3 + 0], xyz[j * 3 + 0]);
            EXPECT_EQ(ref[j * 3 + 1], xyz[j * 3 + 1]);
            EXPECT_EQ(ref[j * 3 + 2], xyz[j * 3 + 2]);
        }
    }

done:
    md_vm_arena_destroy(arena);
    md_file_close(file);
    xdrfile_close(xdr);
}

#if FULL_TEST
UTEST(xtc, big) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));

    const str_t path = STR_LIT("/home/robin/data/PROD_r2.part0001.xtc");
    md_file_o*  file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);

    XDRFILE* xdr = xdrfile_open(path.ptr, "r");

    ASSERT_TRUE(file);
    ASSERT_TRUE(xdr);

    md_array(int64_t) frame_offsets = 0;
    md_array(double)  frame_times = 0;

    const size_t num_atoms  = 5502934;
    const size_t num_frames = 241;

    const size_t coord_size = num_atoms * 3 * sizeof(float);
    float *ref = (float*)md_vm_arena_push(arena, coord_size);
    float *xyz = (float*)md_vm_arena_push(arena, coord_size);

    md_xtc_read_frame_offsets_and_times(file, &frame_offsets, &frame_times, arena);
    size_t xtc_num_frames = frame_offsets ? md_array_size(frame_offsets) - 1 : 0;
    EXPECT_EQ(xtc_num_frames, num_frames);

    int natoms, step;
    float time, box[3][3];

    for (size_t i = 0; i < 1; ++i) {
        md_file_seek(file, frame_offsets[i], MD_FILE_BEG);
        EXPECT_TRUE(md_xtc_read_frame_header(file, &natoms, &step, &time, box));
        EXPECT_TRUE(md_xtc_read_frame_coords(file, (float*)xyz, num_atoms));

        static const size_t xtc_header_size = 52;
        xdr_seek(xdr, frame_offsets[i] + xtc_header_size, SEEK_SET);
        int ncoord = natoms;
        float prec;
        if (xdrfile_decompress_coord_float(ref, &ncoord, &prec, xdr) != natoms) {
            MD_LOG_ERROR("Error reading coordinates from XDR file\n");
            goto done;
        }
        EXPECT_EQ(ncoord, num_atoms);

        for (int j = 0; j < num_atoms; ++j) {
            //EXPECT_EQ(ref[j * 3 + 0], xyz[j * 3 + 0]);
            //EXPECT_EQ(ref[j * 3 + 1], xyz[j * 3 + 1]);
            //EXPECT_EQ(ref[j * 3 + 2], xyz[j * 3 + 2]);
        }
    }
    
done:
    md_vm_arena_destroy(arena);
    md_file_close(file);
    xdrfile_close(xdr);
}



UTEST(xtc, H1N1) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));

    const str_t path = STR_LIT("E:/data/md/H1N1/H1N1-Mich2015-TRAJECTORY-not_water_not_ions-sk100.xtc");
    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);

    XDRFILE* xdr = xdrfile_open(path.ptr, "r");

    ASSERT_TRUE(file);
    ASSERT_TRUE(xdr);

    md_array(int64_t) frame_offsets = 0;
    md_array(double)  frame_times = 0;

    const size_t num_atoms  = 14009213;
    const size_t num_frames = 71;

    const size_t coord_size = num_atoms * 3 * sizeof(float);
    float *ref = (float*)md_vm_arena_push(arena, coord_size);
    float *xyz = (float*)md_vm_arena_push(arena, coord_size);

    md_xtc_read_frame_offsets_and_times(file, &frame_offsets, &frame_times, arena);
    size_t xtc_num_frames = frame_offsets ? md_array_size(frame_offsets) - 1 : 0;
    EXPECT_EQ(xtc_num_frames, num_frames);

    int natoms, step;
    float time, box[3][3];

    for (size_t i = 0; i < 1; ++i) {
        md_file_seek(file, frame_offsets[i], MD_FILE_BEG);
        EXPECT_TRUE(md_xtc_read_frame_header(file, &natoms, &step, &time, box));
        EXPECT_TRUE(md_xtc_read_frame_coords(file, xyz, num_atoms));

        static const size_t xtc_header_size = 52;
        xdr_seek(xdr, frame_offsets[i] + xtc_header_size, SEEK_SET);
        int ncoord = natoms;
        float prec;
        if (xdrfile_decompress_coord_float(ref, &ncoord, &prec, xdr) != natoms) {
            MD_LOG_ERROR("Error reading coordinates from XDR file\n");
            goto done;
        }
        EXPECT_EQ(ncoord, num_atoms);

        for (int j = 0; j < num_atoms; ++j) {
            EXPECT_EQ(ref[j * 3 + 0], xyz[j * 3 + 0]);
            EXPECT_EQ(ref[j * 3 + 1], xyz[j * 3 + 1]);
            EXPECT_EQ(ref[j * 3 + 2], xyz[j * 3 + 2]);
        }
    }

done:
    md_vm_arena_destroy(arena);
    md_file_close(file);
    xdrfile_close(xdr);
}
#endif

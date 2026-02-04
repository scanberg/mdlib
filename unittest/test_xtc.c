#include "utest.h"
#include <string.h>

#include <md_xtc.h>
#include <md_trajectory.h>
#include <md_system.h>
#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_log.h>

#include <xdrfile_xtc.h>

#define FULL_TEST 1
#define ONLY_ONE_FRAME 1
#define XTC_COORD_OFFSET 52

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

static const int num_bits[] = {64, 48, 7, 1, 2, 64, 5, 32, 8, 55, 8, 55, 55, 48, 8, 4, 1, 1, 2, 3, 7, 64, 64, 64, 64, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32 ,33, 35, 41, 44, 51, 55, 59, 63, 63, 63, 1, 2, 3,4};

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

    for (size_t i = 0; i < ARRAY_SIZE(num_bits); ++i) {
        uint64_t ref = decodebits(buf, num_bits[i]);
        uint64_t val = br_read(&r, num_bits[i]);
        EXPECT_EQ(val, ref);
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

static inline uint64_t extract_bits_57_be(
    const uint64_t* base,
    size_t bit_offset,
    size_t bit_length)
{
    size_t wi = bit_offset >> 6;
    size_t sh = bit_offset & 63;
    
    uint64_t word = base[wi];
    
    // Left shift to align, then right shift to extract
    return (word << sh) >> (64 - bit_length);
}

static inline uint32_t extract_bits_be_32(
    const uint64_t* base,
    size_t bit_offset,
    size_t bit_length)   // 1..32
{
    size_t wi = bit_offset >> 6;
    size_t sh = bit_offset & 63;
    
    uint64_t word = base[wi];
    
    // Extract bit_length bits starting at bit position sh
    return (uint32_t)((word << sh) >> (64 - bit_length));
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
    memcpy(&raw, &base[byte_offset], sizeof(uint64_t));
    
    // Byteswap to host endian (assuming host is little-endian)
    uint64_t swapped = BSWAP64(raw);
    
    // Extract the bits (they're now at the top after bswap)
    return (swapped << bit_in_byte) >> (64 - bit_length);
}

static inline uint64_t extract_bits_be_raw_64(
    const uint8_t* base,
    size_t bit_offset,
    size_t bit_length)   // 1..64
{
    size_t byte_offset = bit_offset >> 3;
    size_t bit_in_byte = bit_offset & 7;
    
    // Always read 16 bytes (two uint64_t)
    uint64_t raw_hi, raw_lo;
    memcpy(&raw_hi, &base[byte_offset], sizeof(uint64_t));
    memcpy(&raw_lo, &base[byte_offset + 8], sizeof(uint64_t));
    
    uint64_t hi = BSWAP64(raw_hi);
    uint64_t lo = BSWAP64(raw_lo);
    
    // Combine (lo part masked out if not needed)
    uint64_t nz_mask = -(uint64_t)(bit_in_byte != 0);
    uint64_t combined = (hi << bit_in_byte) | ((lo >> (64 - bit_in_byte)) & nz_mask);
    
    return combined >> (64 - bit_length);
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
    memcpy(&raw, &base[byte_offset], sizeof(uint64_t));
    
    uint64_t swapped = BSWAP64(raw);
    return (uint32_t)((swapped << bit_in_byte) >> (64 - bit_length));
}

static inline uint64_t extract_bits_be(
    const uint64_t* base,
    size_t bit_offset,
    size_t bit_length)   // 1..64
{
    size_t wi = bit_offset >> 6;
    size_t sh = bit_offset & 63;
    
    uint64_t hi = base[wi];
    uint64_t lo = base[wi + 1];
    
    // Alternative: use conditional move (may compile to cmov)
    uint64_t lo_part = (sh == 0) ? 0 : (lo >> (64 - sh));
    uint64_t win_hi = (hi << sh) | lo_part;
    
    return win_hi >> (64 - bit_length);
}

uint64_t extract_packed_le(
    const uint64_t* base,
    size_t bit_offset,
    size_t nbits)
{
    size_t wi  = bit_offset >> 6;
    size_t sh  = bit_offset & 63;

    uint64_t hi = BSWAP64(base[wi]);
    uint64_t lo = BSWAP64(base[wi + 1]);

    // Build 128-bit MSB-first window (hi:lo)
    uint64_t rshift  = (64 - sh) & 63;
    uint64_t lo_part = lo >> rshift;
    uint64_t nz_mask = -(uint64_t)(sh != 0);

    uint64_t win_hi = (hi << sh) | (lo_part & nz_mask);

    // Field is MSB-aligned in win_hi, move to LSB
    uint64_t x = win_hi >> (64 - nbits);

    // Byte packing
    uint64_t byte_count  = (nbits + 7) >> 3;
    uint64_t pad_bits    = (byte_count << 3) - nbits;
    uint64_t final_shift = 64 - (byte_count << 3);

    x >>= pad_bits;

    return BSWAP64(x) >> final_shift;
}

UTEST(xtc, decode_bits) {
    size_t temp_pos = md_temp_get_pos();
    size_t size_in_bytes = sizeof(int) * 1024;
    int* buf = md_temp_push(3 * sizeof(int) + size_in_bytes);
    buf[0] = buf[1] = buf[2] = 0;
    srand(0);
    for (int i = 3; i < 1024; ++i) {
        buf[i] = rand() << 16 | rand();
    }

    const uint64_t* buf_raw = (const uint64_t*)&buf[3];

    uint64_t* buf_be = md_temp_push(size_in_bytes);
    size_t num_u64 = size_in_bytes / sizeof(uint64_t);
    // Preform pre-BE swaps
    for (size_t i = 0; i < num_u64; ++i) {
        buf_be[i] = BSWAP64(buf_raw[i]);
    }

    br_t r = {0};
    br_init(&r, buf_raw, 512 * sizeof(int));

    brn_t brn;
	brn_init(&brn, buf_raw, 512 * sizeof(int));
    
    size_t bit_offset = 0;

    for (size_t i = 0; i < ARRAY_SIZE(num_bits); ++i) {
        int num_of_bits = num_bits[i];
           
        int fullbytes = num_of_bits >> 3;
        int partbits  = num_of_bits &  7;

        uint64_t v = 0;
        int i = 0;
        for (; i < fullbytes; i++) {
            uint64_t ibyte = decodebits(buf, 8);
            v |= ibyte << (8 * i);
        }

        if (partbits) {
            v |= ((uint64_t) decodebits(buf, partbits)) << (8 * i);
        }

        int part_bits = num_of_bits & 7;
        uint64_t big_shift = 64 - ((num_of_bits + 7) & ~7);  // Align to the next multiple of 8
        uint64_t sml_shift = (8 - partbits) & 7;             // Avoid branch, ensures zero when partbits == 0
        uint64_t part_mask = part_bits ? (1ull << part_bits) - 1 : 0xFF;

        uint64_t w = br_read(&r, num_of_bits);
        // Remove unnecessary masking by directly applying shifts
        uint64_t t = (w << sml_shift) & ~0xFFull | (w & part_mask);
        // Byte-swap and shift to align the result correctly
        uint64_t q = BSWAP64(t) >> big_shift;

        uint64_t x = extract_bits_be(buf_be, bit_offset, num_of_bits);
        uint64_t y = extract_bits_be_raw_64((const uint8_t*)buf_raw, bit_offset, num_of_bits);

		uint64_t z = brn_read_u64(&brn, num_of_bits);

        EXPECT_EQ(x, w);
        EXPECT_EQ(y, w);
        EXPECT_EQ(v, q);
        if (v != q) {
            printf("Error with N=%i\n", num_of_bits);
        }

        bit_offset += num_of_bits;
    }

    md_temp_set_pos_back(temp_pos);
}

UTEST(xtc, trajectory_i) {
    const str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.xtc");
    md_trajectory_i* traj = md_xtc_trajectory_create(path, md_get_heap_allocator(), MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE);
    ASSERT_TRUE(traj);

    const size_t num_atoms  = md_trajectory_num_atoms(traj);
    const size_t num_frames = md_trajectory_num_frames(traj);

    EXPECT_EQ(num_atoms, 1336);
    EXPECT_EQ(num_frames, 501);

    const size_t mem_size = num_atoms * 3 * sizeof(float);
    void* mem_ptr = md_alloc(md_get_heap_allocator(), mem_size);
    float *x = (float*)mem_ptr;
    float *y = (float*)mem_ptr + num_atoms * 1;
    float *z = (float*)mem_ptr + num_atoms * 2;

    md_trajectory_frame_header_t header;

	size_t N = ONLY_ONE_FRAME ? 1 : num_frames;

    for (size_t i = 0; i < N; ++i) {
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

        xdr_seek(xdr, frame_offsets[i] + XTC_COORD_OFFSET, SEEK_SET);
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

    const str_t path = STR_LIT("E:/data/md/big/PROD_r2.part0001.xtc");
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

	size_t N = ONLY_ONE_FRAME ? 1 : num_frames;

    for (size_t i = 0; i < N; ++i) {
        md_file_seek(file, frame_offsets[i], MD_FILE_BEG);
        EXPECT_TRUE(md_xtc_read_frame_header(file, &natoms, &step, &time, box));
        EXPECT_TRUE(md_xtc_read_frame_coords(file, (float*)xyz, num_atoms));

        xdr_seek(xdr, frame_offsets[i] + XTC_COORD_OFFSET, SEEK_SET);
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

	size_t N = ONLY_ONE_FRAME ? 1 : num_frames;

    for (size_t i = 0; i < N; ++i) {
        md_file_seek(file, frame_offsets[i], MD_FILE_BEG);
        EXPECT_TRUE(md_xtc_read_frame_header(file, &natoms, &step, &time, box));
        EXPECT_TRUE(md_xtc_read_frame_coords(file, xyz, num_atoms));

        xdr_seek(xdr, frame_offsets[i] + XTC_COORD_OFFSET, SEEK_SET);
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

	size_t N = ONLY_ONE_FRAME ? 1 : num_frames;

    for (size_t i = 0; i < N; ++i) {
        md_file_seek(file, frame_offsets[i], MD_FILE_BEG);
        EXPECT_TRUE(md_xtc_read_frame_header(file, &natoms, &step, &time, box));
        EXPECT_TRUE(md_xtc_read_frame_coords(file, xyz, num_atoms));

        xdr_seek(xdr, frame_offsets[i] + XTC_COORD_OFFSET, SEEK_SET);
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

UTEST(xtc, f1) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));

    const str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/xtc/f1.xtc");
    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);

    XDRFILE* xdr = xdrfile_open(path.ptr, "r");

    ASSERT_TRUE(file);
    ASSERT_TRUE(xdr);

    md_array(int64_t) frame_offsets = 0;
    md_array(double)  frame_times = 0;

    const size_t num_atoms  = 29;
    const size_t num_frames = 1;

    const size_t coord_size = num_atoms * 3 * sizeof(float);
    float *ref = (float*)md_vm_arena_push(arena, coord_size);
    float *xyz = (float*)md_vm_arena_push(arena, coord_size);

    md_xtc_read_frame_offsets_and_times(file, &frame_offsets, &frame_times, arena);
    size_t xtc_num_frames = frame_offsets ? md_array_size(frame_offsets) - 1 : 0;
    EXPECT_EQ(xtc_num_frames, num_frames);

    int natoms, step;
    float time, box[3][3];

	size_t N = ONLY_ONE_FRAME ? 1 : num_frames;

    for (size_t i = 0; i < N; ++i) {
        md_file_seek(file, frame_offsets[i], MD_FILE_BEG);
        EXPECT_TRUE(md_xtc_read_frame_header(file, &natoms, &step, &time, box));
        EXPECT_TRUE(md_xtc_read_frame_coords(file, xyz, num_atoms));

        xdr_seek(xdr, frame_offsets[i] + XTC_COORD_OFFSET, SEEK_SET);
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

UTEST(xtc, pep) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));

    const str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/xtc/pep.xtc");
    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);

    XDRFILE* xdr = xdrfile_open(path.ptr, "r");

    ASSERT_TRUE(file);
    ASSERT_TRUE(xdr);

    md_array(int64_t) frame_offsets = 0;
    md_array(double)  frame_times = 0;

    const size_t num_atoms  = 29;
    const size_t num_frames = 21;

    const size_t coord_size = num_atoms * 3 * sizeof(float);
    float *ref = (float*)md_vm_arena_push(arena, coord_size);
    float *xyz = (float*)md_vm_arena_push(arena, coord_size);

    md_xtc_read_frame_offsets_and_times(file, &frame_offsets, &frame_times, arena);
    size_t xtc_num_frames = frame_offsets ? md_array_size(frame_offsets) - 1 : 0;
    EXPECT_EQ(xtc_num_frames, num_frames);

    int natoms, step;
    float time, box[3][3];

	size_t N = ONLY_ONE_FRAME ? 1 : num_frames;

    for (size_t i = 0; i < N; ++i) {
        md_file_seek(file, frame_offsets[i], MD_FILE_BEG);
        EXPECT_TRUE(md_xtc_read_frame_header(file, &natoms, &step, &time, box));
        EXPECT_TRUE(md_xtc_read_frame_coords(file, xyz, num_atoms));

        xdr_seek(xdr, frame_offsets[i] + XTC_COORD_OFFSET, SEEK_SET);
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

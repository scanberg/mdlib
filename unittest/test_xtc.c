#include "utest.h"
#include <string.h>
#include <stdio.h>

#include <md_xtc.h>
#include <md_gro.h>
#include <md_script.h>
#include <md_system.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_log.h>

#include <xdrfile_xtc.h>

#include "run_check.h"

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

// The stream is the file's bytes, which need not sit on an 8 byte boundary
static inline uint64_t load_qword(const uint64_t* p) {
    uint64_t v;
    memcpy(&v, p, sizeof(v));
    return v;
}

static void br_init(br_t* r, const uint64_t* stream, size_t num_qwords) {
    ASSERT(num_qwords >= 2);
    r->stream = stream + 2;
    r->data   = load_qword(stream + 0);
    r->next   = load_qword(stream + 1);
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
        uint64_t data = load_qword(r->stream++);
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

    uint64_t data = load_qword(r->stream++);
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
        uint64_t data = load_qword(r->stream++);
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
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    int* buf = md_temp_alloc(temp, sizeof(int) * (1024 + 3));
    buf[0] = buf[1] = buf[2] = 0;
    srand(0);
    for (int i = 3; i < 1024; ++i) {
        buf[i] = (int)((unsigned)rand() << 16 | (unsigned)rand());
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

    md_temp_end(temp);
}

UTEST(xtc, decode_bits) {
    md_temp_scope_t temp = md_temp_begin();
    int* buf = md_temp_alloc(temp, sizeof(int) * (1024 + 3));
    buf[0] = buf[1] = buf[2] = 0;
    srand(0);
    for (int i = 3; i < 1024; ++i) {
        buf[i] = (int)((unsigned)rand() << 16 | (unsigned)rand());
    }

    br_t r = {0};
    br_init(&r, (const uint64_t*)&buf[3], 512 * sizeof(int));

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

        uint64_t big_shift = 64 - ((num_of_bits + 7) & ~7);  // Align to the next multiple of 8
        uint64_t sml_shift = (8 - partbits) & 7;             // Avoid branch, ensures zero when partbits == 0
        uint64_t part_mask = partbits ? (1ull << partbits) - 1 : 0xFF;

        uint64_t next = ALIGN_TO(num_of_bits, 8);
        uint64_t u = br_peek(&r, num_of_bits);
        uint64_t l = u << sml_shift;
        uint64_t j = BSWAP64(u);
        uint64_t m = BSWAP64(l);
        uint64_t k = j >> big_shift;

        uint64_t w = br_read(&r, num_of_bits);
        // Remove unnecessary masking by directly applying shifts
        uint64_t t = (w << sml_shift) & ~0xFFull | (w & part_mask);
        // Byte-swap and shift to align the result correctly
        uint64_t q = BSWAP64(t) >> big_shift;

        EXPECT_EQ(v, q);
        if (v != q) {
            printf("Error with N=%i\n", num_of_bits);
        }
    }

    md_temp_end(temp);
}

// Every frame of the file extracts, with no structure to go with it: the run is the file alone.
UTEST(xtc, run_every_frame) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    const str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/catalyst.xtc");
    ASSERT_TRUE(md_xtc_system_publish_run(&sys, path, STR_LIT("run/c"), MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    EXPECT_EQ(501u, run_num_frames(&sys, STR_LIT("run/c")));
    EXPECT_EQ(1336u, run_num_atoms(&sys, STR_LIT("run/c")));

    md_system_state_t st = {.alloc = arena};
    md_system_state_init(&st, 1336);
    const str_t paths[] = { STR_INIT("atom/position"), STR_INIT("unitcell") };
    md_system_extract_t* ex = md_system_extract_begin(&sys, STR_LIT("run/c"), paths, 2, md_get_heap_allocator());
    ASSERT_TRUE(ex != NULL);
    for (int64_t i = 0; i < 501; ++i) {
        EXPECT_TRUE(md_system_extract_frame(ex, i, &st));
    }
    md_system_extract_end(ex);
    md_vm_arena_destroy(arena);
}

UTEST(xtc, catalyst) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp);

    const str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/catalyst.xtc");
    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));

    XDRFILE* xdr = xdrfile_open(path.ptr, "r");

    ASSERT_TRUE(xdr);

    md_array(int64_t) frame_offsets = 0;
    md_array(double)  frame_times = 0;

    const size_t num_atoms  = 1336;
    const size_t num_frames = 501;

    const size_t coord_size = num_atoms * 3 * sizeof(float);
    float *ref = (float*)md_temp_alloc(temp, coord_size);
    float *xyz = (float*)md_temp_alloc(temp, coord_size);

    md_xtc_read_frame_offsets_and_times(file, &frame_offsets, &frame_times, arena);
    size_t xtc_num_frames = frame_offsets ? md_array_size(frame_offsets) - 1 : 0;
    EXPECT_EQ(xtc_num_frames, num_frames);

    md_array(uint8_t) frame_data = 0;
    md_xtc_header_t xtc_header = { 0 };

    for (size_t i = 0; i < num_frames; ++i) {
        md_file_offset_t frame_beg = frame_offsets[i];
        md_file_offset_t frame_end = frame_offsets[i + 1];
        size_t frame_size = frame_end - frame_beg;
        md_file_seek(file, frame_beg, MD_FILE_BEG);

        md_array_ensure(frame_data, frame_size + MD_XTC_DECODE_PADDING, arena);
        size_t read_bytes = md_file_read_at(file, frame_beg, frame_data, frame_size);
        EXPECT_TRUE(md_xtc_decode_frame_data(frame_data, read_bytes, &xtc_header, xyz, num_atoms));

        static const size_t xtc_header_size = 52;
        xdr_seek(xdr, frame_beg + xtc_header_size, SEEK_SET);
        int ncoord = xtc_header.natoms;
        float prec;
        if (xdrfile_decompress_coord_float(ref, &ncoord, &prec, xdr) != xtc_header.natoms) {
            MD_LOG_ERROR("Error reading coordinates from XDR file\n");
            goto done;
        }
        ASSERT_EQ(ncoord, num_atoms);

        for (int j = 0; j < num_atoms; ++j) {
            ASSERT_EQ(ref[j * 3 + 0], xyz[j * 3 + 0]);
            ASSERT_EQ(ref[j * 3 + 1], xyz[j * 3 + 1]);
            ASSERT_EQ(ref[j * 3 + 2], xyz[j * 3 + 2]);
        }
    }
    
done:
    md_temp_end(temp);
    md_file_close(&file);
    xdrfile_close(xdr);
}

#if FULL_TEST
UTEST(xtc, big) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp_scope);

    const str_t path = STR_INIT("E:/data/md/big/PROD_r2.part0001.xtc");
    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));

    XDRFILE* xdr = xdrfile_open(path.ptr, "r");

    ASSERT_TRUE(xdr);

    md_array(int64_t) frame_offsets = 0;
    md_array(double)  frame_times = 0;

    const size_t num_atoms  = 5502934;
    const size_t num_frames = 241;

    const size_t coord_size = num_atoms * 3 * sizeof(float);
    float *ref = (float*)md_temp_alloc(temp_scope, coord_size);
    float *xyz = (float*)md_temp_alloc(temp_scope, coord_size);

    md_xtc_read_frame_offsets_and_times(file, &frame_offsets, &frame_times, arena);
    size_t xtc_num_frames = frame_offsets ? md_array_size(frame_offsets) - 1 : 0;
    EXPECT_EQ(xtc_num_frames, num_frames);

    md_array(uint8_t) frame_data = 0;
    md_xtc_header_t xtc_header = { 0 };

    for (size_t i = 0; i < 1; ++i) {
        md_file_offset_t frame_beg = frame_offsets[i];
        md_file_offset_t frame_end = frame_offsets[i + 1];
        size_t frame_size = frame_end - frame_beg;
        md_file_seek(file, frame_beg, MD_FILE_BEG);

        md_array_ensure(frame_data, frame_size + MD_XTC_DECODE_PADDING, arena);
        size_t read_bytes = md_file_read_at(file, frame_beg, frame_data, frame_size);
        EXPECT_TRUE(md_xtc_decode_frame_data(frame_data, read_bytes, &xtc_header, xyz, num_atoms));

        static const size_t xtc_header_size = 52;
        xdr_seek(xdr, frame_offsets[i] + xtc_header_size, SEEK_SET);
        int ncoord = xtc_header.natoms;
        float prec;
        if (xdrfile_decompress_coord_float(ref, &ncoord, &prec, xdr) != xtc_header.natoms) {
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
    md_temp_end(temp_scope);
    md_file_close(&file);
    xdrfile_close(xdr);
}

UTEST(xtc, amyloid) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp_scope);

    const str_t path = STR_INIT("E:/data/md/amyloid-6T/prod-centered.xtc");
    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));

    XDRFILE* xdr = xdrfile_open(path.ptr, "r");

    ASSERT_TRUE(xdr);

    md_array(int64_t) frame_offsets = 0;
    md_array(double)  frame_times = 0;

    const size_t num_atoms  = 161271;
    const size_t num_frames = 5701;

    const size_t coord_size = num_atoms * 3 * sizeof(float);
    float *ref = (float*)md_temp_alloc(temp_scope, coord_size);
    float *xyz = (float*)md_temp_alloc(temp_scope, coord_size);

    md_xtc_read_frame_offsets_and_times(file, &frame_offsets, &frame_times, arena);
    size_t xtc_num_frames = frame_offsets ? md_array_size(frame_offsets) - 1 : 0;
    EXPECT_EQ(xtc_num_frames, num_frames);

    md_array(uint8_t) frame_data = 0;
    md_xtc_header_t xtc_header = { 0 };

    for (size_t i = 0; i < 5; ++i) {
        md_file_offset_t frame_beg = frame_offsets[i];
        md_file_offset_t frame_end = frame_offsets[i + 1];
        size_t frame_size = frame_end - frame_beg;
        md_file_seek(file, frame_beg, MD_FILE_BEG);

        static const size_t xtc_header_size = 52;
        xdr_seek(xdr, frame_offsets[i] + xtc_header_size, SEEK_SET);
        int ncoord = xtc_header.natoms;
        float prec;
        if (xdrfile_decompress_coord_float(ref, &ncoord, &prec, xdr) != xtc_header.natoms) {
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
    md_temp_end(temp_scope);
    md_file_close(&file);
    xdrfile_close(xdr);
}

UTEST(xtc, H1N1) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp);

    const str_t path = STR_INIT("E:/data/md/H1N1/H1N1-Mich2015-TRAJECTORY-not_water_not_ions-sk100.xtc");
    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));

    XDRFILE* xdr = xdrfile_open(path.ptr, "r");

    ASSERT_TRUE(xdr);

    md_array(int64_t) frame_offsets = 0;
    md_array(double)  frame_times = 0;

    const size_t num_atoms  = 14009213;
    const size_t num_frames = 71;

    const size_t coord_size = num_atoms * 3 * sizeof(float);
    float *ref = (float*)md_temp_alloc(temp, coord_size);
    float *xyz = (float*)md_temp_alloc(temp, coord_size);

    md_xtc_read_frame_offsets_and_times(file, &frame_offsets, &frame_times, arena);
    size_t xtc_num_frames = frame_offsets ? md_array_size(frame_offsets) - 1 : 0;
    EXPECT_EQ(xtc_num_frames, num_frames);

    md_array(uint8_t) frame_data = 0;
    md_xtc_header_t xtc_header = { 0 };

    for (size_t i = 0; i < 1; ++i) {
        md_file_offset_t frame_beg = frame_offsets[i];
        md_file_offset_t frame_end = frame_offsets[i + 1];
        size_t frame_size = frame_end - frame_beg;
        md_file_seek(file, frame_beg, MD_FILE_BEG);

        static const size_t xtc_header_size = 52;
        xdr_seek(xdr, frame_offsets[i] + xtc_header_size, SEEK_SET);
        int ncoord = xtc_header.natoms;
        float prec;
        if (xdrfile_decompress_coord_float(ref, &ncoord, &prec, xdr) != xtc_header.natoms) {
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
    md_temp_end(temp);
    md_file_close(&file);
    xdrfile_close(xdr);
}
#endif

// ### RUN ###

#define XTC_RUN STR_LIT("run/catalyst")

static const str_t xtc_coord_paths[] = { STR_INIT("atom/position"), STR_INIT("unitcell") };

// Recorded from the trajectory reader this replaced
static const run_ref_t xtc_refs[] = {
    { 0,   {43228.999, 40956.939, 69523.6185},   {24.6800003, 36.5499992, 20.0299988}, {36.8499985, 23.8400002, 72.2799988}, {61.8499985, 66.1999969, 103.310005, 0, 0, 0} },
    { 250, {41309.5191, 46800.459, 72117.2884},  {24.789999, 36.6100006, 20.0200005},  {20.5799999, 36.9899979, 77.1199951}, {61.8499985, 66.1999969, 103.310005, 0, 0, 0} },
    { 500, {42644.869, 42075.7291, 70371.9886},  {24.7199993, 36.5499992, 20.1000004}, {39.579998, 38.3400002, 74.7200012},  {61.8499985, 66.1999969, 103.310005, 0, 0, 0} },
};

static bool xtc_load_catalyst(md_system_t* sys, md_allocator_i* arena) {
    const str_t xtc = STR_INIT(MD_UNITTEST_DATA_DIR "/catalyst.xtc");
    sys->alloc = arena;
    md_system_state_t sys_state = {.alloc = arena};
    return md_gro_system_init_from_file(sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.gro"), NULL) &&
        md_xtc_system_publish_run(sys, xtc, XTC_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE);
}

// The run holds the file's frames as attributes, and a state extracted from it holds what the
// trajectory reader it replaced gave for the same frame.
UTEST(xtc, run_matches_reference) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {0};
    ASSERT_TRUE(xtc_load_catalyst(&sys, arena));
    run_check_refs(utest_result, &sys, XTC_RUN, 501, 1336, xtc_refs, ARRAY_SIZE(xtc_refs));

    const md_attributes_t* t = &sys.attributes;
    const size_t F = 501;
    const size_t N = sys.atom.count;

    const md_attribute_t* time = md_attributes_find(t, STR_LIT("run/catalyst/time"));
    const md_attribute_t* step = md_attributes_find(t, STR_LIT("run/catalyst/step"));
    const md_attribute_t* cell = md_attributes_find(t, STR_LIT("run/catalyst/unitcell"));
    const md_attribute_t* pos  = md_attributes_find(t, STR_LIT("run/catalyst/atom/position"));
    ASSERT_TRUE(time && step && cell && pos);
    EXPECT_EQ(F, step->format.shape[0]);
    EXPECT_EQ(3u, pos->format.components);
    EXPECT_EQ(time, md_attributes_axis(t, pos));
    EXPECT_EQ(time, md_attributes_axis(t, cell));
    EXPECT_TRUE(md_unit_equal(time->unit, md_unit_picosecond()));
    for (size_t i = 1; i < F; ++i) {
        EXPECT_LT(((const double*)time->data)[i - 1], ((const double*)time->data)[i]);
        EXPECT_LT(((const int64_t*)step->data)[i - 1], ((const int64_t*)step->data)[i]);
    }

    md_system_state_t got = {.alloc = arena};
    md_system_state_init(&got, N);

    md_system_extract_t* ex = md_system_extract_begin(&sys, XTC_RUN, xtc_coord_paths, ARRAY_SIZE(xtc_coord_paths), md_get_heap_allocator());
    ASSERT_TRUE(ex != NULL);
    ASSERT_TRUE(md_system_extract_frame(ex, 7, &got));
    EXPECT_EQ(7.0, got.frame);
    // Out of range.
    EXPECT_FALSE(md_system_extract_frame(ex, (int64_t)F, &got));
    md_system_extract_end(ex);

    // One atom of one frame, through the attribute and without a context.
    float xyz[3];
    md_attribute_slice_t one = md_attribute_slice_2((uint32_t)F / 2, 7);
    ASSERT_TRUE(run_extract_one(&got, &sys, XTC_RUN, (int64_t)F / 2));
    ASSERT_EQ(md_attribute_extract_f32(xyz, 3, pos, one, md_unit_none()), 3u);
    EXPECT_EQ(got.xyz[7].x, xyz[0]);
    EXPECT_EQ(got.xyz[7].y, xyz[1]);
    EXPECT_EQ(got.xyz[7].z, xyz[2]);

    // And in nanometer, converted rather than reinterpreted.
    ASSERT_EQ(md_attribute_extract_f32(xyz, 3, pos, one, md_unit_nanometer()), 3u);
    EXPECT_NEAR(got.xyz[7].x * 0.1f, xyz[0], 1.0e-5f);

    // Every frame at once is exactly what the virtual attribute exists to avoid.
    float* all = md_alloc(arena, F * N * 3 * sizeof(float));
    EXPECT_EQ(md_attribute_extract_f32(all, F * N * 3, pos, md_attribute_slice_all(), md_unit_none()), 0u);

    // The cell and the frame alone: only the cell asked for, into a state without coordinates.
    const str_t cell_only[] = { STR_INIT("unitcell") };
    ex = md_system_extract_begin(&sys, XTC_RUN, cell_only, 1, md_get_heap_allocator());
    ASSERT_TRUE(ex != NULL);
    md_system_state_t meta = {0};
    ASSERT_TRUE(md_system_extract_frame(ex, 3, &meta));
    EXPECT_EQ(3.0, meta.frame);
    EXPECT_NE(0u, meta.unitcell.flags);
    md_system_extract_end(ex);

    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

UTEST(xtc, run_extracts_other_attributes_into_the_state) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {0};
    ASSERT_TRUE(xtc_load_catalyst(&sys, arena));
    const size_t F = run_num_frames(&sys, XTC_RUN);
    const double* frame_times = (const double*)md_attributes_find(&sys.attributes, STR_LIT("run/catalyst/time"))->data;

    // Twice the rate of the trajectory, as an energy file often is.
    const uint32_t R = (uint32_t)(2 * F - 1);
    double*  obs_time  = md_alloc(arena, R * sizeof(double));
    double*  obs_value = md_alloc(arena, R * sizeof(double));
    int32_t* obs_label = md_alloc(arena, F * 3 * sizeof(int32_t));
    for (uint32_t r = 0; r < R; ++r) {
        obs_time[r]  = (r % 2 == 0) ? frame_times[r / 2] : 0.5 * (frame_times[r / 2] + frame_times[r / 2 + 1]);
        obs_value[r] = 2.0 * obs_time[r];
    }
    for (size_t i = 0; i < F * 3; ++i) obs_label[i] = (int32_t)i;

    md_attributes_t* t = &sys.attributes;
    ASSERT_NE(md_attributes_create(t, &(md_attribute_desc_t){ .path = STR_INIT("run/catalyst/obs/time"),
        .format = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = {R}}, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .unit = md_unit_picosecond(), .data = obs_time, .byte_size = R * sizeof(double)}), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(md_attributes_create(t, &(md_attribute_desc_t){ .path = STR_INIT("run/catalyst/obs/value"),
        .format = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = {R}}, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .unit = md_unit_kelvin(), .data = obs_value, .byte_size = R * sizeof(double)}), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(md_attributes_create(t, &(md_attribute_desc_t){ .path = STR_INIT("run/catalyst/label"),
        .format = {.type = MD_ATTRIBUTE_TYPE_I32, .components = 1, .rank = 2, .shape = {(uint32_t)F, 3}}, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .unit = md_unit_none(), .data = obs_label, .byte_size = F * 3 * sizeof(int32_t)}), MD_ATTRIBUTE_INVALID);

    const str_t paths[] = { STR_INIT("obs/value"), STR_INIT("label") };
    md_system_extract_t* ex = md_system_extract_begin(&sys, XTC_RUN, paths, ARRAY_SIZE(paths), md_get_heap_allocator());
    ASSERT_TRUE(ex != NULL);

    md_system_state_t st = {.alloc = arena};
    md_system_state_init(&st, 0);
    for (int64_t f = 0; f < (int64_t)F; f += 11) {
        ASSERT_TRUE(md_system_extract_frame(ex, f, &st));
        const md_attribute_t* v = md_attributes_find(&st.attributes, STR_LIT("obs/value"));
        const md_attribute_t* l = md_attributes_find(&st.attributes, STR_LIT("label"));
        ASSERT_TRUE(v && l);
        // The value at the frame's own time, not at row f of the finer axis.
        EXPECT_EQ(0u, v->format.rank);
        EXPECT_EQ(0u, v->flags & MD_ATTRIBUTE_FLAG_TEMPORAL);
        EXPECT_TRUE(md_unit_equal(v->unit, md_unit_kelvin()));
        EXPECT_NEAR(2.0 * frame_times[f], ((const double*)v->data)[0], 1.0e-9);
        // Integers stay integers.
        EXPECT_EQ(MD_ATTRIBUTE_TYPE_I32, l->format.type);
        EXPECT_EQ(1u, l->format.rank);
        EXPECT_EQ(3u, l->format.shape[0]);
        EXPECT_EQ((int32_t)(f * 3 + 2), ((const int32_t*)l->data)[2]);
    }

    // Removing what the context reads fails the next extract rather than reading freed storage.
    md_attributes_remove(t, md_attributes_find(t, STR_LIT("run/catalyst/label"))->id);
    EXPECT_FALSE(md_system_extract_frame(ex, 0, &st));
    md_system_extract_end(ex);

    md_system_state_free(&st);
    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

// A script evaluated along the run gives the same values in one range as split over ranges the way
// a pool splits it, each range with its own extraction context.
UTEST(xtc, script_evaluates_along_the_run) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {0};
    ASSERT_TRUE(xtc_load_catalyst(&sys, arena));
    const uint32_t F = (uint32_t)run_num_frames(&sys, XTC_RUN);

    md_script_ir_t* ir = md_script_ir_create(arena);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("d = distance(1, 200);"), &sys, NULL));

    md_script_eval_t* a = md_script_eval_create(F, ir, arena);
    md_script_eval_t* b = md_script_eval_create(F, ir, arena);
    ASSERT_TRUE(md_script_eval_frame_range(a, ir, &sys, XTC_RUN, 0, F));
    ASSERT_TRUE(md_script_eval_frame_range(b, ir, &sys, XTC_RUN, 0, F / 3));
    ASSERT_TRUE(md_script_eval_frame_range(b, ir, &sys, XTC_RUN, F / 3, F));

    const md_attribute_t* da = md_attributes_find(md_script_eval_attributes(a), STR_LIT("script/d"));
    const md_attribute_t* db = md_attributes_find(md_script_eval_attributes(b), STR_LIT("script/d"));
    ASSERT_TRUE(da && db);
    EXPECT_EQ(0, MEMCMP(da->data, db->data, md_attribute_byte_size(&da->format)));

    md_script_eval_free(a);
    md_script_eval_free(b);
    md_script_ir_free(ir);
    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

// A trajectory that does not fit the system is refused, rather than published and read out of step.
UTEST(xtc, run_refuses_a_different_system) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t sys_state = {.alloc = arena};
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/water.gro"), NULL));
    EXPECT_FALSE(md_xtc_system_publish_run(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.xtc"), XTC_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

// A run is its positions: a frame axis alone gives nothing to extract, and neither does no run.
UTEST(xtc, run_needs_positions) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    sys.attributes.alloc = arena;
    const double times[3] = { 0.0, 1.0, 2.0 };
    ASSERT_NE(md_attributes_create(&sys.attributes, &(md_attribute_desc_t){
        .path = STR_INIT("run/bare/time"),
        .format = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = {3}},
        .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_picosecond(),
        .data = times, .byte_size = sizeof(times)}), MD_ATTRIBUTE_INVALID);

    EXPECT_TRUE(md_system_extract_begin(&sys, STR_LIT("run/bare"), xtc_coord_paths, 2, md_get_heap_allocator()) == NULL);
    EXPECT_TRUE(md_system_extract_begin(&sys, (str_t){0}, xtc_coord_paths, 2, md_get_heap_allocator()) == NULL);
    EXPECT_TRUE(md_system_extract_begin(&sys, STR_LIT("run/none"), xtc_coord_paths, 2, md_get_heap_allocator()) == NULL);

    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

#if MD_PLATFORM_UNIX
// The point of the context: the file is opened once and kept. Removing it from the directory after
// the first frame leaves the open file readable on unix, so every later frame still extracts - which
// it could not if the file were opened again per frame, as it is without a context.
UTEST(xtc, run_context_keeps_the_file_open) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    const str_t src = STR_INIT(MD_UNITTEST_DATA_DIR "/catalyst.xtc");
    const str_t tmp = STR_INIT("/tmp/md_unittest_run_context.xtc");

    // A private copy to remove.
    md_file_t in = {0}, out = {0};
    ASSERT_TRUE(md_file_open(&in, src, MD_FILE_READ));
    ASSERT_TRUE(md_file_open(&out, tmp, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
    const size_t size = md_file_size(in);
    void* bytes = md_alloc(arena, size);
    ASSERT_EQ(md_file_read(in, bytes, size), size);
    ASSERT_EQ(md_file_write(out, bytes, size), size);
    md_file_close(&in);
    md_file_close(&out);

    md_system_t sys = {.alloc = arena};
    md_system_state_t sys_state = {.alloc = arena};
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.gro"), NULL));
    ASSERT_TRUE(md_xtc_system_publish_run(&sys, tmp, XTC_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));

    md_system_state_t st = {.alloc = arena};
    md_system_state_init(&st, sys.atom.count);
    md_system_extract_t* ex = md_system_extract_begin(&sys, XTC_RUN, xtc_coord_paths, 2, md_get_heap_allocator());
    ASSERT_TRUE(ex != NULL);
    ASSERT_TRUE(md_system_extract_frame(ex, 0, &st));

    ASSERT_EQ(0, remove(tmp.ptr));
    EXPECT_TRUE(md_system_extract_frame(ex, 1, &st));
    EXPECT_TRUE(md_system_extract_frame(ex, 400, &st));
    md_system_extract_end(ex);

    // Without the context each frame opens the file, which is gone.
    const md_attribute_t* pos = md_attributes_find(&sys.attributes, STR_LIT("run/catalyst/atom/position"));
    float* xyz = md_alloc(arena, sys.atom.count * 3 * sizeof(float));
    md_attribute_slice_t s1 = md_attribute_slice_1(1);
    EXPECT_EQ(0u, md_attribute_extract_f32(xyz, sys.atom.count * 3, pos, s1, md_unit_none()));

    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}
#endif

UTEST(xtc, atom_count_from_the_first_frame) {
    EXPECT_EQ(1336u, md_xtc_file_atom_count(STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.xtc")));
    EXPECT_EQ(0u, md_xtc_file_atom_count(STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.gro")));
}

// The two ends of the coordinate compression which ordinary files do not reach, as files written here from the bytes
// below. Each has 12 atoms, since a frame of 9 or fewer is stored uncompressed.
//
// xtc_planar (written by MDAnalysis 2): two frames, every atom at the same y and z, so the range of those is a single
// value. The three way division of the small integers then divides by 1, which libdivide's branchfree divider
// refuses with an abort: a planar system (a 2D material, a membrane model on a grid) brought the whole program down.
static const uint8_t xtc_planar[248] = {
    0x00, 0x00, 0x07, 0xcb, 0x00, 0x00, 0x00, 0x0c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x40, 0xa0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x40, 0xa0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x40, 0xa0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0c, 0x44, 0x7a, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x07, 0xd0, 0x00, 0x00, 0x05, 0xdc, 0x00, 0x00, 0x04, 0x4c, 0x00, 0x00, 0x07, 0xd0,
    0x00, 0x00, 0x05, 0xdc, 0x00, 0x00, 0x00, 0x14, 0x00, 0x00, 0x00, 0x1e, 0x00, 0x10, 0xb2, 0x08,
    0xb2, 0x04, 0x45, 0x86, 0x2f, 0x43, 0x2b, 0xaf, 0x88, 0x0b, 0xc5, 0x2c, 0x04, 0x00, 0xe4, 0x23,
    0x94, 0x5d, 0x35, 0x84, 0x99, 0x25, 0x79, 0x53, 0x04, 0x40, 0x00, 0x00, 0x00, 0x00, 0x07, 0xcb,
    0x00, 0x00, 0x00, 0x0c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0xa0, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0xa0, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0xa0, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x0c, 0x44, 0x7a, 0x00, 0x00, 0x00, 0x00, 0x00, 0x32, 0x00, 0x00, 0x07, 0xd0,
    0x00, 0x00, 0x05, 0xdc, 0x00, 0x00, 0x04, 0x7e, 0x00, 0x00, 0x07, 0xd0, 0x00, 0x00, 0x05, 0xdc,
    0x00, 0x00, 0x00, 0x14, 0x00, 0x00, 0x00, 0x1e, 0x00, 0x10, 0xb2, 0x08, 0xb2, 0x04, 0x45, 0x86,
    0x2f, 0x43, 0x2b, 0xaf, 0x88, 0x0b, 0xc5, 0x2c, 0x04, 0x00, 0xe4, 0x23, 0x94, 0x5d, 0x35, 0x84,
    0x99, 0x25, 0x79, 0x53, 0x04, 0x40, 0x00, 0x00,
};

// xtc_wide_range (written by hand, and read back by a port of the xdrfile reference): x spans 300 nm at precision
// 1e5, a range above 0xffffff, for which the frame states each coordinate as three plain bit fields rather than one
// packed integer. Those were read with the byte reordering of the packed form, which put an atom 2450 A from where
// it was.
static const uint8_t xtc_wide_range[192] = {
    0x00, 0x00, 0x07, 0xcb, 0x00, 0x00, 0x00, 0x0c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x40, 0xa0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x40, 0xa0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x40, 0xa0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0c, 0x47, 0xc3, 0x50, 0x00, 0xff, 0x1b, 0x1e, 0x40,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x86, 0xa0, 0x00, 0xbe, 0xbc, 0x20, 0x00, 0x04, 0x32, 0x38,
    0x00, 0x09, 0xeb, 0x10, 0x00, 0x00, 0x00, 0x09, 0x00, 0x00, 0x00, 0x62, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x09, 0x89, 0x68, 0x03, 0x0d, 0x40, 0x61, 0xa8, 0x09, 0x89, 0x68, 0x03,
    0x0d, 0x40, 0x61, 0xa8, 0x07, 0x27, 0x0e, 0x02, 0x49, 0xf0, 0x49, 0x3e, 0x04, 0xc4, 0xb4, 0x01,
    0x86, 0xa0, 0x30, 0xd4, 0x02, 0xfa, 0xf0, 0x80, 0xf4, 0x24, 0x1e, 0x84, 0x81, 0xc9, 0xc3, 0x80,
    0x92, 0x7c, 0x12, 0x4f, 0x81, 0x0b, 0x07, 0x60, 0x55, 0x73, 0x0a, 0xae, 0x60, 0x98, 0x96, 0x80,
    0x30, 0xd4, 0x06, 0x1a, 0x80, 0x55, 0xd4, 0xa8, 0x1b, 0x77, 0x43, 0x6e, 0xe8, 0x2f, 0xaf, 0x08,
    0x0f, 0x42, 0x41, 0xe8, 0x48, 0x1a, 0x39, 0xde, 0x08, 0x64, 0x71, 0x0c, 0x8e, 0x00, 0x00, 0x00,
};

static bool xtc_write_bytes(str_t path, const uint8_t* bytes, size_t size) {
    md_file_t f = {0};
    if (!md_file_open(&f, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) return false;
    const bool ok = md_file_write(f, bytes, size) == size;
    md_file_close(&f);
    return ok;
}

UTEST(xtc, planar_frames) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    const str_t path = STR_LIT("md_unittest_planar.xtc");
    ASSERT_TRUE(xtc_write_bytes(path, xtc_planar, sizeof(xtc_planar)));
    ASSERT_TRUE(md_xtc_system_publish_run(&sys, path, STR_LIT("run/p"), MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    ASSERT_EQ(2u, run_num_frames(&sys, STR_LIT("run/p")));
    md_system_state_t st = {.alloc = arena};
    md_system_state_init(&st, 12);
    for (int f = 0; f < 2; ++f) {
        ASSERT_TRUE(run_extract_one(&st, &sys, STR_LIT("run/p"), f));
        for (int i = 0; i < 12; ++i) {
            EXPECT_NEAR(i + 0.5f * f, st.xyz[i].x, 1.0e-3f);
            EXPECT_NEAR(20.0f, st.xyz[i].y, 1.0e-3f);
            EXPECT_NEAR(15.0f, st.xyz[i].z, 1.0e-3f);
        }
    }
    remove(path.ptr);
    md_vm_arena_destroy(arena);
}

UTEST(xtc, ranges_beyond_24_bits) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    const str_t path = STR_LIT("md_unittest_wide_range.xtc");
    ASSERT_TRUE(xtc_write_bytes(path, xtc_wide_range, sizeof(xtc_wide_range)));
    ASSERT_TRUE(md_xtc_system_publish_run(&sys, path, STR_LIT("run/w"), MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    md_system_state_t st = {.alloc = arena};
    md_system_state_init(&st, 12);
    ASSERT_TRUE(run_extract_one(&st, &sys, STR_LIT("run/w"), 0));
    for (int i = 0; i < 12; ++i) {
        // nm: (-150 + 25 i, 0.25 i, 1 + 0.5 i)
        EXPECT_NEAR(-1500.0f + 250.0f * i, st.xyz[i].x, 1.0e-3f);
        EXPECT_NEAR(2.5f * i,              st.xyz[i].y, 1.0e-3f);
        EXPECT_NEAR(10.0f + 5.0f * i,      st.xyz[i].z, 1.0e-3f);
    }
    remove(path.ptr);
    md_vm_arena_destroy(arena);
}

// The frames above restated as GROMACS 2023 writes a system of more than 298261617 atoms (and its own tests any
// system): magic 2023, and the byte count of the compressed coordinates as a 64 bit integer, high half first.
// Nothing else differs, so the same coordinates must come out.
static size_t xtc_to_magic_2023(uint8_t* dst, const uint8_t* src, size_t size) {
    size_t i = 0, o = 0;
    while (i + 92 <= size) {
        memcpy(dst + o, src + i, 88);
        dst[o + 0] = 0; dst[o + 1] = 0; dst[o + 2] = 0x07; dst[o + 3] = 0xe7;   // 2023
        const uint32_t n = ((uint32_t)src[i + 88] << 24) | ((uint32_t)src[i + 89] << 16) | ((uint32_t)src[i + 90] << 8) | src[i + 91];
        memset(dst + o + 88, 0, 4);
        memcpy(dst + o + 92, src + i + 88, 4);
        const size_t padded = (n + 3) & ~(size_t)3;
        memcpy(dst + o + 96, src + i + 92, padded);
        i += 92 + padded;
        o += 96 + padded;
    }
    return o;
}

UTEST(xtc, magic_2023) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    const str_t path = STR_LIT("md_unittest_magic_2023.xtc");
    uint8_t bytes[512];
    struct { const uint8_t* src; size_t size; int frames; } files[2] = {
        { xtc_planar,     sizeof(xtc_planar),     2 },
        { xtc_wide_range, sizeof(xtc_wide_range), 1 },
    };
    for (int k = 0; k < 2; ++k) {
        const size_t size = xtc_to_magic_2023(bytes, files[k].src, files[k].size);
        ASSERT_EQ(files[k].size + 4 * (size_t)files[k].frames, size);
        ASSERT_TRUE(xtc_write_bytes(path, bytes, size));
        EXPECT_EQ(12u, md_xtc_file_atom_count(path));
        md_system_t sys = {.alloc = arena};
        ASSERT_TRUE(md_xtc_system_publish_run(&sys, path, STR_LIT("run/n"), MD_RUN_FLAG_DISABLE_CACHE_WRITE));
        ASSERT_EQ((size_t)files[k].frames, run_num_frames(&sys, STR_LIT("run/n")));
        md_system_state_t st = {.alloc = arena};
        md_system_state_init(&st, 12);
        for (int f = 0; f < files[k].frames; ++f) {
            ASSERT_TRUE(run_extract_one(&st, &sys, STR_LIT("run/n"), f));
            for (int i = 0; i < 12; ++i) {
                if (k == 0) {
                    EXPECT_NEAR(i + 0.5f * f, st.xyz[i].x, 1.0e-3f);
                    EXPECT_NEAR(15.0f, st.xyz[i].z, 1.0e-3f);
                } else {
                    EXPECT_NEAR(-1500.0f + 250.0f * i, st.xyz[i].x, 1.0e-3f);
                    EXPECT_NEAR(10.0f + 5.0f * i, st.xyz[i].z, 1.0e-3f);
                }
            }
        }
        md_system_free(&sys);
    }
    remove(path.ptr);
    md_vm_arena_destroy(arena);
}

#include "utest.h"
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <hdf5.h>

#include <md_hdf5.h>
#include <md_hdf5_direct.h>
#include <md_attributes.h>
#include "../src/md_attributes_internal.h"
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_os.h>
#include <core/md_str.h>

// Datasets stored every way md_hdf5_direct decodes, and some it leaves to HDF5, each read back frame by frame and in
// pieces, and held against what HDF5 itself reads.

#define HD_PATH "md_unittest_hdf5_direct.h5"

typedef struct hd_case_t {
    const char* name;
    int      rank;
    hsize_t  dims[4];
    hsize_t  chunk[4];          // all zero for contiguous
    hid_t    file_type;         // a predefined type, copied
    bool     shuffle;
    int      deflate;
    bool     fletcher;
    bool     fletcher_first;    // the checksum before the other filters rather than after
    bool     nbit;              // a filter not decoded here
    bool     leave_unwritten;   // chunks of the last frame never written
    bool     unfiltered_chunk;  // the first chunk written past its filters
    bool     shuffled_only;     // ... past its deflate only: shuffled, not deflated
    bool     direct;            // what the storage is expected to be read as
} hd_case_t;

static double hd_value(size_t flat) {
    // Fits every type below: small integers, exact in float
    return (double)((flat * 7) % 251);
}

static size_t hd_count(const hsize_t* dims, int rank) {
    size_t n = 1;
    for (int i = 0; i < rank; ++i) n *= (size_t)dims[i];
    return n;
}

static bool hd_write(const hd_case_t* c) {
    const hid_t file = H5Fcreate(HD_PATH, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) return false;
    const size_t n = hd_count(c->dims, c->rank);
    double* values = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; ++i) values[i] = hd_value(i);

    hsize_t max[4];
    for (int i = 0; i < c->rank; ++i) max[i] = c->dims[i];
    const bool chunked = c->chunk[0] != 0;
    if (chunked) max[0] = H5S_UNLIMITED;
    const hid_t space = H5Screate_simple(c->rank, c->dims, max);
    const hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    if (chunked) {
        H5Pset_chunk(dcpl, c->rank, c->chunk);
        if (c->fletcher && c->fletcher_first) H5Pset_fletcher32(dcpl);
        if (c->nbit)     H5Pset_nbit(dcpl);
        if (c->shuffle)  H5Pset_shuffle(dcpl);
        if (c->deflate)  H5Pset_deflate(dcpl, (unsigned)c->deflate);
        if (c->fletcher && !c->fletcher_first) H5Pset_fletcher32(dcpl);
    }
    H5Pset_fill_value(dcpl, H5T_NATIVE_DOUBLE, &(double){ 99.0 });
    const hid_t type = H5Tcopy(c->file_type);
    const hid_t d = H5Dcreate2(file, "data", type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    bool ok = d >= 0;
    if (ok && c->leave_unwritten) {
        // All frames but the last
        hsize_t start[4] = {0}, count[4];
        for (int i = 0; i < c->rank; ++i) count[i] = c->dims[i];
        count[0] -= 1;
        const hid_t fs = H5Dget_space(d);
        const hid_t ms = H5Screate_simple(c->rank, count, NULL);
        H5Sselect_hyperslab(fs, H5S_SELECT_SET, start, NULL, count, NULL);
        ok = H5Dwrite(d, H5T_NATIVE_DOUBLE, ms, fs, H5P_DEFAULT, values) >= 0;
        H5Sclose(ms);
        H5Sclose(fs);
    } else if (ok) {
        ok = H5Dwrite(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values) >= 0;
    }
    if (ok && c->unfiltered_chunk) {
        // The first chunk again, as it is, with every filter marked as skipped
        const size_t elem = H5Tget_size(type);
        const size_t chunk_elems = hd_count(c->chunk, c->rank);
        uint8_t* raw = calloc(chunk_elems, elem);
        hsize_t stride[4] = {0};
        hsize_t s = 1;
        for (int i = c->rank - 1; i >= 0; --i) { stride[i] = s; s *= c->dims[i]; }
        for (size_t e = 0; e < chunk_elems; ++e) {
            // The element's place in the dataset
            size_t rem = e, flat = 0;
            bool inside = true;
            for (int i = c->rank - 1; i >= 0; --i) {
                const size_t k = rem % c->chunk[i];
                rem /= c->chunk[i];
                inside &= k < c->dims[i];
                flat += k * stride[i];
            }
            if (!inside) continue;
            double tmp = hd_value(flat);   // converted in place: a double holds every type here
            H5Tconvert(H5T_NATIVE_DOUBLE, type, 1, &tmp, NULL, H5P_DEFAULT);
            MD_MEMCPY(raw + e * elem, &tmp, elem);
        }
        const hsize_t offset[4] = {0};
        uint32_t mask = 0xffffffffu;
        if (c->shuffled_only) {
            // Shuffle is filter 0 and deflate filter 1: the bytes shuffled, deflate marked as skipped
            uint8_t* shuffled = malloc(chunk_elems * elem);
            for (size_t b = 0; b < elem; ++b) for (size_t e = 0; e < chunk_elems; ++e) shuffled[b * chunk_elems + e] = raw[e * elem + b];
            MD_MEMCPY(raw, shuffled, chunk_elems * elem);
            free(shuffled);
            mask = 1u << 1;
        }
        ok = H5Dwrite_chunk(d, H5P_DEFAULT, mask, offset, chunk_elems * elem, raw) >= 0;
        free(raw);
    }
    if (d >= 0) H5Dclose(d);
    H5Tclose(type);
    H5Pclose(dcpl);
    H5Sclose(space);
    H5Fclose(file);
    free(values);
    return ok;
}

// What HDF5 reads, all of it, as doubles
static double* hd_reference(const hd_case_t* c) {
    const size_t n = hd_count(c->dims, c->rank);
    double* out = malloc(n * sizeof(double));
    const hid_t file = H5Fopen(HD_PATH, H5F_ACC_RDONLY, H5P_DEFAULT);
    const hid_t d = H5Dopen2(file, "data", H5P_DEFAULT);
    H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out);
    H5Dclose(d);
    H5Fclose(file);
    return out;
}

UTEST(hdf5_direct, storage_of_every_kind) {
    const hd_case_t cases[] = {
        { "contiguous f32",           3, {5, 13, 3},    {0},          H5T_IEEE_F32LE, .direct = true },
        { "contiguous f64 be",        2, {4, 9},        {0},          H5T_IEEE_F64BE, .direct = true },
        { "frame chunks, deflate",    3, {5, 13, 3},    {1, 13, 3},   H5T_IEEE_F32LE, .shuffle = true, .deflate = 4, .direct = true },
        { "3 frames a chunk",         3, {7, 13, 3},    {3, 13, 3},   H5T_IEEE_F32LE, .shuffle = true, .deflate = 1, .direct = true },
        { "split atoms and xyz",      3, {5, 13, 3},    {2, 5, 2},    H5T_IEEE_F32LE, .shuffle = true, .deflate = 6, .direct = true },
        { "rank 4",                   4, {3, 4, 5, 3},  {2, 3, 2, 3}, H5T_IEEE_F64LE, .deflate = 2, .direct = true },
        { "rank 1",                   1, {11},          {4},          H5T_IEEE_F64LE, .deflate = 2, .direct = true },
        { "u16 be, fletcher32",       3, {5, 13, 3},    {1, 13, 3},   H5T_STD_U16BE,  .shuffle = true, .deflate = 3, .fletcher = true, .direct = true },
        { "i32, fletcher32 first",    2, {6, 17},       {2, 8},       H5T_STD_I32LE,  .fletcher = true, .fletcher_first = true, .deflate = 1, .direct = true },
        { "i8, raw chunks",           2, {6, 17},       {4, 5},       H5T_STD_I8LE,   .direct = true },
        { "a chunk unfiltered",       3, {5, 13, 3},    {1, 13, 3},   H5T_IEEE_F32LE, .shuffle = true, .deflate = 4, .unfiltered_chunk = true, .direct = true },
        { "a chunk not deflated",     3, {5, 13, 3},    {1, 13, 3},   H5T_IEEE_F32LE, .shuffle = true, .deflate = 4, .unfiltered_chunk = true, .shuffled_only = true, .direct = true },
        { "a frame never written",    3, {5, 13, 3},    {1, 13, 3},   H5T_IEEE_F32LE, .deflate = 4, .leave_unwritten = true, .direct = true },
        { "nbit: through HDF5",       3, {5, 13, 3},    {1, 13, 3},   H5T_STD_I32LE,  .nbit = true, .deflate = 1, .direct = false },
    };
    for (size_t ci = 0; ci < MD_COUNTOF(cases); ++ci) {
        const hd_case_t* c = &cases[ci];
        ASSERT_TRUE(hd_write(c));
        double* ref = hd_reference(c);

        const md_hdf5_lock_t lock = md_hdf5_lock();
        const hid_t file = H5Fopen(HD_PATH, H5F_ACC_RDONLY, H5P_DEFAULT);
        const hid_t dset = H5Dopen2(file, "data", H5P_DEFAULT);
        md_hdf5_direct_t* d = md_hdf5_direct_describe(file, dset, md_get_heap_allocator());
        H5Dclose(dset);
        H5Fclose(file);
        md_hdf5_unlock(lock);
        ASSERT_TRUE(d != NULL);
        EXPECT_EQ(c->direct, md_hdf5_direct_is_direct(d));
        if (c->direct != md_hdf5_direct_is_direct(d)) printf("%s: mode %u\n", c->name, d->mode);

        const size_t F = (size_t)c->dims[0];
        const size_t row = hd_count(c->dims + 1, c->rank - 1);
        float* out = malloc(row * sizeof(float));
        int wrong = 0;
        for (int pass = 0; pass < 2; ++pass) {
            // Without an io, and with one that keeps the decoded chunks
            md_attribute_io_t io = { .alloc = md_get_heap_allocator() };
            md_attribute_io_t* pio = pass ? &io : NULL;
            for (size_t f = F; f-- > 0;) {
                // The whole row, then pieces of it of every length from every start
                if (!md_hdf5_direct_read_f32(out, d, STR_LIT(HD_PATH), "data", f, 0, row, 2.0, pio)) { wrong += 1000; continue; }
                for (size_t i = 0; i < row; ++i) wrong += out[i] != (float)(2.0 * ref[f * row + i]);
                for (size_t first = 0; first < row; first += 3) {
                    for (size_t count = 1; first + count <= row; count += 4) {
                        if (!md_hdf5_direct_read_f32(out, d, STR_LIT(HD_PATH), "data", f, first, count, 1.0, pio)) { wrong += 1000; continue; }
                        for (size_t i = 0; i < count; ++i) wrong += out[i] != (float)ref[f * row + first + i];
                    }
                }
            }
            md_attribute_io_close_all(&io);
        }
        EXPECT_EQ(0, wrong);
        if (wrong) printf("%s: %d wrong\n", c->name, wrong);
        // Out of range
        EXPECT_FALSE(md_hdf5_direct_read_f32(out, d, STR_LIT(HD_PATH), "data", F, 0, 1, 1.0, NULL));
        EXPECT_FALSE(md_hdf5_direct_read_f32(out, d, STR_LIT(HD_PATH), "data", 0, row, 1, 1.0, NULL));

        // Frame spans: contiguous and whole frame chunks have them
        int64_t offsets[16], sizes[16];
        md_hdf5_direct_frame_spans(offsets, sizes, d, F);
        const bool whole_frames = c->chunk[0] == 0 || (c->rank == 1 || hd_count(c->chunk + 1, c->rank - 1) >= row);
        if (c->direct && whole_frames && !c->leave_unwritten) {
            for (size_t f = 0; f < F; ++f) EXPECT_GT(offsets[f], 0);
        }
        free(out);
        free(ref);
        md_free(md_get_heap_allocator(), d, (size_t)d->size);
    }
    remove(HD_PATH);
}

UTEST(hdf5_direct, block_cache) {
    md_attribute_io_t io = { .alloc = md_get_heap_allocator() };
    bool hit = true;
    const uint64_t a = md_attribute_io_block_key(STR_LIT("a.h5"), 100);
    const uint64_t b = md_attribute_io_block_key(STR_LIT("a.h5"), 200);
    EXPECT_NE(a, b);
    EXPECT_NE(a, md_attribute_io_block_key(STR_LIT("b.h5"), 100));
    uint8_t* p = md_attribute_io_block(&io, a, 64, &hit);
    ASSERT_TRUE(p != NULL);
    EXPECT_FALSE(hit);
    p[0] = 42;
    EXPECT_TRUE(md_attribute_io_block(&io, a, 64, &hit) == p);
    EXPECT_TRUE(hit);
    // Another size under the key is another block
    md_attribute_io_block(&io, a, 32, &hit);
    EXPECT_FALSE(hit);
    // Kept while fewer than the slots are in use, the least recently used given up after
    for (int i = 0; i < 2 * MD_ATTRIBUTE_IO_MAX_BLOCKS; ++i) md_attribute_io_block(&io, md_attribute_io_block_key(STR_LIT("c.h5"), i), 16, &hit);
    md_attribute_io_block(&io, b, 16, &hit);
    EXPECT_FALSE(hit);
    EXPECT_TRUE(md_attribute_io_block(&io, b, 16, &hit) != NULL);
    EXPECT_TRUE(hit);
    md_attribute_io_block_forget(&io, b);
    md_attribute_io_block(&io, b, 16, &hit);
    EXPECT_FALSE(hit);
    EXPECT_TRUE(md_attribute_io_block(NULL, b, 16, &hit) == NULL);
    md_attribute_io_close_all(&io);
}

// HDF5 1.14.6's H5_checksum_fletcher32, as it is, to hold md_hdf5_fletcher32 against
static uint32_t hd_fletcher_ref(const void* _data, size_t _len) {
    const uint8_t* data = (const uint8_t*)_data;
    size_t len = _len / 2;
    uint32_t sum1 = 0, sum2 = 0;
    while (len) {
        size_t tlen = len > 360 ? 360 : len;
        len -= tlen;
        do {
            sum1 += (uint32_t)(((uint16_t)data[0]) << 8) | ((uint16_t)data[1]);
            data += 2;
            sum2 += sum1;
        } while (--tlen);
        sum1 = (sum1 & 0xffff) + (sum1 >> 16);
        sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    }
    if (_len % 2) {
        sum1 += (uint32_t)(((uint16_t)*data) << 8);
        sum2 += sum1;
        sum1 = (sum1 & 0xffff) + (sum1 >> 16);
        sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    }
    sum1 = (sum1 & 0xffff) + (sum1 >> 16);
    sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    return (sum2 << 16) | sum1;
}

UTEST(hdf5_direct, fletcher32_as_hdf5) {
    const size_t cap = 3 * 4096 * 2 + 77;
    uint8_t* d = malloc(cap);
    uint32_t seed = 12345;
    int wrong = 0, checks = 0;
    // Lengths around the block of 4096 words and the fold of 360, random and odd ones; random bytes, all 0xff (a sum
    // that is a multiple of 65535), all zero, and mostly zero
    size_t lengths[64];
    size_t num_lengths = 0;
    const size_t edges[] = { 1, 2, 3, 31, 32, 33, 719, 720, 721, 8191, 8192, 8193, 8194, 16384, 16385, 24576, 24577 };
    for (size_t i = 0; i < MD_COUNTOF(edges); ++i) lengths[num_lengths++] = edges[i];
    while (num_lengths < MD_COUNTOF(lengths)) { seed = seed * 1664525u + 1013904223u; lengths[num_lengths++] = 1 + seed % cap; }
    for (int fill = 0; fill < 4; ++fill) {
        for (size_t i = 0; i < cap; ++i) {
            seed = seed * 1664525u + 1013904223u;
            d[i] = fill == 0 ? (uint8_t)(seed >> 24) : fill == 1 ? 0xff : fill == 2 ? 0 : ((seed >> 24) % 7 == 0 ? (uint8_t)(seed >> 16) : 0);
        }
        for (size_t l = 0; l < num_lengths; ++l) {
            // From an offset that is not aligned, too
            for (size_t off = 0; off < 2; ++off) {
                const size_t len = MD_MIN(lengths[l], cap - off);
                wrong += md_hdf5_fletcher32(d + off, len) != hd_fletcher_ref(d + off, len);
                checks += 1;
            }
        }
    }
    EXPECT_EQ(0, wrong);
    EXPECT_GT(checks, 500);
    free(d);
}

// A damaged chunk is refused, as HDF5 refuses it; a checksum as HDF5 before 1.6.3 wrote it is accepted
UTEST(hdf5_direct, fletcher32_damage) {
    for (int fletcher_first = 0; fletcher_first < 2; ++fletcher_first) {
        // The checksum last, of the deflated bytes as stored (h5py's order); or first, of the data, which is then
        // shuffled (netCDF-4's order, without the deflate: a damaged deflate stream is zlib's to find)
        const hd_case_t c = fletcher_first
            ? (hd_case_t){ "checksum first", 3, {4, 600, 3}, {1, 600, 3}, H5T_IEEE_F32LE, .shuffle = true, .fletcher = true, .fletcher_first = true, .direct = true }
            : (hd_case_t){ "checksum last",  3, {4, 600, 3}, {1, 600, 3}, H5T_IEEE_F32LE, .shuffle = true, .deflate = 1, .fletcher = true, .direct = true };
        ASSERT_TRUE(hd_write(&c));

        md_hdf5_lock_t lock = md_hdf5_lock();
        hid_t file = H5Fopen(HD_PATH, H5F_ACC_RDONLY, H5P_DEFAULT);
        hid_t dset = H5Dopen2(file, "data", H5P_DEFAULT);
        md_hdf5_direct_t* d = md_hdf5_direct_describe(file, dset, md_get_heap_allocator());
        H5Dclose(dset);
        H5Fclose(file);
        md_hdf5_unlock(lock);
        ASSERT_TRUE(d && md_hdf5_direct_is_direct(d));

        const size_t row = 600 * 3;
        float* out = malloc(row * sizeof(float));
        ASSERT_TRUE(md_hdf5_direct_read_f32(out, d, STR_LIT(HD_PATH), "data", 2, 0, row, 1.0, NULL));

        // The stored checksum of frame 2, and a byte of its chunk
        const md_hdf5_direct_chunk_t e = d->table[2];
        uint8_t* chunk = malloc(e.size);
        md_file_t f = {0};
        ASSERT_TRUE(md_file_open(&f, STR_LIT(HD_PATH), MD_FILE_READ | MD_FILE_WRITE));
        ASSERT_EQ((size_t)e.size, md_file_read_at(f, e.offset, chunk, e.size));

        // The pre 1.6.3 checksum, each half's bytes swapped: where the checksum is the last four bytes stored
        if (!fletcher_first) {
            const uint8_t* sum = chunk + e.size - 4;
            const uint8_t swapped[4] = { sum[1], sum[0], sum[3], sum[2] };
            md_file_write_at(f, e.offset + e.size - 4, swapped, 4);
        }
        md_file_close(&f);
        EXPECT_TRUE(md_hdf5_direct_read_f32(out, d, STR_LIT(HD_PATH), "data", 2, 0, row, 1.0, NULL));

        // One byte in the middle of the chunk flipped: refused, by HDF5 as well
        ASSERT_TRUE(md_file_open(&f, STR_LIT(HD_PATH), MD_FILE_READ | MD_FILE_WRITE));
        uint8_t b = chunk[e.size / 2] ^ 0x10;
        md_file_write_at(f, e.offset + e.size / 2, &b, 1);
        md_file_close(&f);
        EXPECT_FALSE(md_hdf5_direct_read_f32(out, d, STR_LIT(HD_PATH), "data", 2, 0, row, 1.0, NULL));
        // The frames of other chunks are unharmed
        EXPECT_TRUE(md_hdf5_direct_read_f32(out, d, STR_LIT(HD_PATH), "data", 1, 0, row, 1.0, NULL));

        lock = md_hdf5_lock();
        file = H5Fopen(HD_PATH, H5F_ACC_RDONLY, H5P_DEFAULT);
        dset = H5Dopen2(file, "data", H5P_DEFAULT);
        float* all = malloc(4 * row * sizeof(float));
        EXPECT_LT(H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, all), 0);
        free(all);
        H5Dclose(dset);
        H5Fclose(file);
        md_hdf5_unlock(lock);

        free(chunk);
        free(out);
        md_free(md_get_heap_allocator(), d, (size_t)d->size);
    }
    remove(HD_PATH);
}

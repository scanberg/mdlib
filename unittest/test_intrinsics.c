#include "utest.h"

#include <core/md_compiler.h>
#include <core/md_common.h>

#if (MD_COMPILER_GCC || MD_COMPILER_CLANG)
#pragma GCC diagnostic ignored "-Wsign-compare"
#endif

UTEST(intrinsics, memset) {
    char data[131];
    MD_MEMSET(data, 4, sizeof(data));

    for (int i = 0; i < MD_COUNTOF(data); ++i) {
        EXPECT_EQ(data[i], 4);
    }
}

UTEST(intrinsics, memcpy) {
    int a[131];
    int b[131];
    for (int i = 0; i < MD_COUNTOF(a); ++i) {
        a[i] = i;
    }
    MD_MEMCPY(b, a, sizeof(b));

    for (int i = 0; i < MD_COUNTOF(a); ++i) {
        EXPECT_EQ(a[i], b[i]);
    }
}

UTEST(intrinsics, test) {
    {
        static const uint32_t mask[5] = {0, 0xFFU, 0xFFFFU, 0xFFFFFFU, 0xFFFFFFFFU};

        EXPECT_EQ(md_clz32(mask[1]), 24);
        EXPECT_EQ(md_clz32(mask[2]), 16);
        EXPECT_EQ(md_clz32(mask[3]), 8);
        EXPECT_EQ(md_clz32(mask[4]), 0);

        EXPECT_EQ(md_popcnt32(mask[0]), 0);
        EXPECT_EQ(md_popcnt32(mask[1]), 8);
        EXPECT_EQ(md_popcnt32(mask[2]), 16);
        EXPECT_EQ(md_popcnt32(mask[3]), 24);
        EXPECT_EQ(md_popcnt32(mask[4]), 32);
    }

    {
        static const uint32_t mask[5] = {0, (1U << 8), (1U << 16), (1U << 24), (1U << 31)};

        EXPECT_EQ(md_bsf32(mask[0]), 0);
        EXPECT_EQ(md_bsf32(mask[1]), 9);
        EXPECT_EQ(md_bsf32(mask[2]), 17);
        EXPECT_EQ(md_bsf32(mask[3]), 25);
        EXPECT_EQ(md_bsf32(mask[4]), 32);

        EXPECT_EQ(md_bsr32(mask[0]), 0);
        EXPECT_EQ(md_bsr32(mask[1]), 9);
        EXPECT_EQ(md_bsr32(mask[2]), 17);
        EXPECT_EQ(md_bsr32(mask[3]), 25);
        EXPECT_EQ(md_bsr32(mask[4]), 32);
    }

    EXPECT_EQ(md_next_power_of_two32(3), 4);
    EXPECT_EQ(md_next_power_of_two32(31), 32);
    EXPECT_EQ(md_next_power_of_two32(32), 32);
    EXPECT_EQ(md_next_power_of_two32(5), 8);
    EXPECT_EQ(md_next_power_of_two32(8000), 8192);

    EXPECT_EQ(md_next_power_of_two64(sizeof(void*)), sizeof(void*));



    uint64_t res = 0;
    res = md_find_first_zero_byte64(65280);
    res = md_find_first_zero_byte64(255);
    res = md_find_first_zero_byte64(65535);


    EXPECT_EQ(md_find_first_zero_byte32(65280),      0); // 1111111100000000
    EXPECT_EQ(md_find_first_zero_byte32(255),        1); // 0000000011111111
    EXPECT_EQ(md_find_first_zero_byte32(65535),      2); // 000000001111111111111111

}

// The bit functions against a loop over the bits, for 0, single bits, all ones and random values: every compiler and
// every CPU gives the same result for every input, 0 included (clz and ctz of 0 are the width)
static uint64_t ref_rng(uint64_t* s) { *s ^= *s << 13; *s ^= *s >> 7; *s ^= *s << 17; return *s; }

static unsigned ref_clz(uint64_t x, unsigned w) { unsigned n = 0; for (int i = (int)w - 1; i >= 0 && !((x >> i) & 1); --i) ++n; return n; }
static unsigned ref_ctz(uint64_t x, unsigned w) { unsigned n = 0; for (unsigned i = 0; i < w && !((x >> i) & 1); ++i) ++n; return n; }
static unsigned ref_pop(uint64_t x) { unsigned n = 0; for (; x; x &= x - 1) ++n; return n; }

static void ref_umul128(uint64_t a, uint64_t b, uint64_t* lo, uint64_t* hi) {
    // Schoolbook on 16 bit digits, sharing nothing with md_umul128: a column holds at most four products below 2^32
    uint64_t r[8] = {0};
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            r[i + j] += ((a >> (16 * i)) & 0xFFFF) * ((b >> (16 * j)) & 0xFFFF);
        }
    }
    uint64_t carry = 0;
    for (int k = 0; k < 8; ++k) {
        r[k] += carry;
        carry = r[k] >> 16;
        r[k] &= 0xFFFF;
    }
    *lo = r[0] | r[1] << 16 | r[2] << 32 | r[3] << 48;
    *hi = r[4] | r[5] << 16 | r[6] << 32 | r[7] << 48;
}

UTEST(intrinsics, bits_defined_everywhere) {
    uint64_t vals[200];
    size_t n = 0;
    vals[n++] = 0;
    vals[n++] = UINT64_MAX;
    for (unsigned i = 0; i < 64; ++i) vals[n++] = (uint64_t)1 << i;
    uint64_t s = 0x9E3779B97F4A7C15ull;
    while (n < MD_COUNTOF(vals)) {
        const uint64_t r = ref_rng(&s);
        vals[n++] = r >> (r & 63);      // all magnitudes
    }

    for (size_t i = 0; i < n; ++i) {
        const uint64_t x = vals[i];
        const uint32_t y = (uint32_t)x;
        EXPECT_EQ(ref_clz(x, 64), md_clz64(x));
        EXPECT_EQ(ref_ctz(x, 64), md_ctz64(x));
        EXPECT_EQ(ref_clz(y, 32), md_clz32(y));
        EXPECT_EQ(ref_ctz(y, 32), md_ctz32(y));
        EXPECT_EQ(ref_pop(x), md_popcnt64(x));
        EXPECT_EQ(ref_pop(y), md_popcnt32(y));
        EXPECT_EQ(x ? ref_ctz(x, 64) + 1 : 0, md_bsf64(x));
        EXPECT_EQ(x ? 64 - ref_clz(x, 64) : 0, md_bsr64(x));
        EXPECT_EQ(y ? ref_ctz(y, 32) + 1 : 0, md_bsf32(y));
        EXPECT_EQ(y ? 32 - ref_clz(y, 32) : 0, md_bsr32(y));

        unsigned zb = 0; while (zb < 8 && ((x >> (8 * zb)) & 0xFF)) ++zb;
        EXPECT_EQ(zb, md_find_first_zero_byte64(x));
        unsigned zn = 0; while (zn < 4 && ((y >> (4 * zn)) & 0xF)) ++zn;
        EXPECT_EQ(zn, md_find_first_zero_nibble16((uint16_t)y));

        uint64_t lo, hi, rlo, rhi;
        const uint64_t b = vals[(i * 7 + 3) % n];
        lo = md_umul128(x, b, &hi);
        ref_umul128(x, b, &rlo, &rhi);
        EXPECT_EQ(rlo, lo);
        EXPECT_EQ(rhi, hi);
    }

    volatile uint64_t c = 40;
    EXPECT_EQ(40u, md_atomic_fetch_add64(&c, 2));
    EXPECT_EQ(42u, md_atomic_load64(&c));
    EXPECT_EQ(42u, md_atomic_fetch_add64(&c, UINT64_MAX));     // minus one
    EXPECT_EQ(41u, md_atomic_load64(&c));
}

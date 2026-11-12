#pragma once

// FIXED WIDTH FIELDS, THREE AT ONCE
//
// md_parse_fixed_f32 (md_parse.h) reads one field of a fixed width column format. The coordinates of the structure
// formats come three at a time, in fields of 8 characters: "%8.3f" for those of pdb and gro, "%8.4f" for the velocities
// of gro. md_parse_fixed_float3_w8 reads three such fields at once (md_simd: SSE/AVX on x86, NEON on arm64), and takes
// exactly what md_parse_fixed_f32 takes, with the same value: where any of the three is laid out otherwise it gives
// false with out untouched, and the caller reads them one by one.
//
// Its own header, so that the includers of md_parse.h do not take md_simd.h with it.

#include <core/md_common.h>
#include <core/md_simd.h>
#include <core/md_str.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

// WHITESPACE, 64 BYTES AT ONCE
//
// md_parse_whitespace_mask64 gives the whitespace of 64 bytes ('\t', '\n', '\r' and ' ', what md_ascii_is_whitespace and
// md_parse_next_token split by), one bit per byte: what a tokenizer takes the tokens of a line from, 64 bytes at a time
// (md_mmcif.c, the rows of _atom_site). p is aligned to 64, so that no load crosses a page: the bytes around a range may
// be read where a page holds bytes of it, which address sanitizers are told not to look at.

#if defined(__clang__) || defined(__GNUC__)
#   define MD_PARSE_SIMD_NO_SANITIZE __attribute__((no_sanitize_address))
#elif defined(_MSC_VER)
#   define MD_PARSE_SIMD_NO_SANITIZE __declspec(no_sanitize_address)
#else
#   define MD_PARSE_SIMD_NO_SANITIZE
#endif

MD_PARSE_SIMD_NO_SANITIZE static MD_FORCE_INLINE uint64_t md_parse_whitespace_mask64(const char* p) {
    MD_ASSERT(((uintptr_t)p & 63) == 0);
    md_128i m[4];
    for (int k = 0; k < 4; ++k) {
        const md_128i v = md_mm_load_si128((const md_128i*)(p + 16 * k));
        m[k] = md_mm_or_si128(md_mm_or_si128(md_mm_cmpeq_epi8(v, md_mm_set1_epi8(' ')),  md_mm_cmpeq_epi8(v, md_mm_set1_epi8('\t'))),
                              md_mm_or_si128(md_mm_cmpeq_epi8(v, md_mm_set1_epi8('\n')), md_mm_cmpeq_epi8(v, md_mm_set1_epi8('\r'))));
    }
    return md_mm_movemask4_epi8(m[0], m[1], m[2], m[3]);
}

// A BYTE, 64 AT ONCE, AND THE SEARCHES BUILT ON IT
//
// md_parse_byte_mask64 gives the bytes equal to c of 64, one bit per byte, as md_parse_whitespace_mask64 does whitespace,
// with p aligned to 64 the same way. On it, in place of memchr (a call, and a setup for any length), the searches of a
// range [beg, end) inline: a count of c, past the n:th, and the lines one after the other. These take the chunks one
// after the other whatever the matches are, rather than searching from where the last match was, which is a chain of
// loads each waiting for the one before: twice as fast over lines of 80 characters. A single search
// (md_parse_find_byte) mostly ends in its first register or two, and takes a register at a time.

MD_PARSE_SIMD_NO_SANITIZE static MD_FORCE_INLINE uint64_t md_parse_byte_mask64(const char* p, char c) {
    MD_ASSERT(((uintptr_t)p & 63) == 0);
#if defined(__AVX2__)
    const md_256i k  = md_mm256_set1_epi8(c);
    const uint64_t lo = (uint32_t)md_mm256_movemask_epi8(md_mm256_cmpeq_epi8(md_mm256_load_si256((const md_256i*)p), k));
    const uint64_t hi = (uint32_t)md_mm256_movemask_epi8(md_mm256_cmpeq_epi8(md_mm256_load_si256((const md_256i*)(p + 32)), k));
    return lo | (hi << 32);
#else
    const md_128i k = md_mm_set1_epi8(c);
    md_128i m[4];
    for (int i = 0; i < 4; ++i) {
        m[i] = md_mm_cmpeq_epi8(md_mm_load_si128((const md_128i*)(p + 16 * i)), k);
    }
    return md_mm_movemask4_epi8(m[0], m[1], m[2], m[3]);
#endif
}

// The chunk of 64 bytes p is in
static MD_FORCE_INLINE const char* md_parse_chunk64(const char* p) {
    return (const char*)((uintptr_t)p & ~(uintptr_t)63);
}

// The bytes equal to c of the register of bytes at p (MD_PARSE_REG_WIDTH of them, p aligned to it), MD_PARSE_REG_BITS bits
// of the mask per byte, the first lowest: what a single search takes, which mostly ends in its first register or two. On
// arm64 the narrowing shift of the comparison (shrn) is what movemask is on x86, four bits per byte.
#if defined(__AVX2__)
#   define MD_PARSE_REG_WIDTH 32
#   define MD_PARSE_REG_BITS  1
#elif defined(MD_SIMD_NEON)
#   define MD_PARSE_REG_WIDTH 16
#   define MD_PARSE_REG_BITS  4
#else
#   define MD_PARSE_REG_WIDTH 16
#   define MD_PARSE_REG_BITS  1
#endif

MD_PARSE_SIMD_NO_SANITIZE static MD_FORCE_INLINE uint64_t md_parse_byte_mask_reg(const char* p, char c) {
#if defined(__AVX2__)
    return (uint32_t)md_mm256_movemask_epi8(md_mm256_cmpeq_epi8(md_mm256_load_si256((const md_256i*)p), md_mm256_set1_epi8(c)));
#elif defined(MD_SIMD_NEON)
    const uint8x16_t eq = vceqq_u8(vld1q_u8((const uint8_t*)p), vdupq_n_u8((uint8_t)c));
    return vget_lane_u64(vreinterpret_u64_u8(vshrn_n_u16(vreinterpretq_u16_u8(eq), 4)), 0);
#else
    return (uint32_t)md_mm_movemask_epi8(md_mm_cmpeq_epi8(md_mm_load_si128((const md_128i*)p), md_mm_set1_epi8(c)));
#endif
}

// The first c of [beg, end), else end: a register at a time
MD_PARSE_SIMD_NO_SANITIZE static inline const char* md_parse_find_byte(const char* beg, const char* end, char c) {
    if (beg >= end) return end;
    const char* a = (const char*)((uintptr_t)beg & ~(uintptr_t)(MD_PARSE_REG_WIDTH - 1));
    uint64_t m = md_parse_byte_mask_reg(a, c) >> ((uintptr_t)(beg - a) * MD_PARSE_REG_BITS);
    if (m) {
        const char* p = beg + md_ctz64(m) / MD_PARSE_REG_BITS;
        return p < end ? p : end;
    }
    for (a += MD_PARSE_REG_WIDTH; a < end; a += MD_PARSE_REG_WIDTH) {
        m = md_parse_byte_mask_reg(a, c);
        if (m) {
            const char* p = a + md_ctz64(m) / MD_PARSE_REG_BITS;
            return p < end ? p : end;
        }
    }
    return end;
}

// The count of c in [beg, end)
MD_PARSE_SIMD_NO_SANITIZE static inline size_t md_parse_count_byte(const char* beg, const char* end, char c) {
    if (beg >= end) return 0;
    const char* chunk = md_parse_chunk64(beg);
    uint64_t m = md_parse_byte_mask64(chunk, c) & (~UINT64_C(0) << (uintptr_t)(beg - chunk));
    size_t n = 0;
    for (;;) {
        const size_t past = (size_t)MD_MAX((ptrdiff_t)0, (chunk + 64) - end);
        if (past) {
            // The last chunk: not past end
            m &= ~UINT64_C(0) >> past;
            return n + (size_t)md_popcnt64(m);
        }
        n += (size_t)md_popcnt64(m);
        chunk += 64;
        if (chunk >= end) return n;   // end at the end of a chunk: none past it is loaded
        m = md_parse_byte_mask64(chunk, c);
    }
}

// Past the n:th c of [beg, end), or end where there are fewer, with *n lessened by those there were
MD_PARSE_SIMD_NO_SANITIZE static inline const char* md_parse_after_nth_byte(const char* beg, const char* end, char c, size_t* n) {
    if (beg >= end || *n == 0) return beg;
    const char* chunk = md_parse_chunk64(beg);
    uint64_t m = md_parse_byte_mask64(chunk, c) & (~UINT64_C(0) << (uintptr_t)(beg - chunk));
    for (;;) {
        const size_t past = (size_t)MD_MAX((ptrdiff_t)0, (chunk + 64) - end);
        if (past) m &= ~UINT64_C(0) >> past;
        const size_t k = (size_t)md_popcnt64(m);
        if (k >= *n) {
            // The n:th of this chunk's: the lowest n - 1 cleared
            for (size_t i = 1; i < *n; ++i) m &= m - 1;
            *n = 0;
            return chunk + md_ctz64(m) + 1;
        }
        *n -= k;
        chunk += 64;
        if (past || chunk >= end) return end;
        m = md_parse_byte_mask64(chunk, c);
    }
}

// The lines of [beg, end) one after the other, each from the mask of the newlines of its chunk alone (its lowest bit,
// cleared): md_parse_lines_next gives the next without its newline (nor a carriage return before it), false after the
// last. The last line of the range need not end in a newline.
typedef struct md_parse_lines_t {
    const char* beg;        // Of the next line
    const char* end;        // Of the range
    const char* chunk;      // The chunk of mask
    uint64_t    mask;       // Its newlines not taken yet
} md_parse_lines_t;

MD_PARSE_SIMD_NO_SANITIZE static inline md_parse_lines_t md_parse_lines(const char* beg, const char* end) {
    md_parse_lines_t it = { .beg = beg, .end = end, .chunk = md_parse_chunk64(beg) };
    it.mask = beg < end ? md_parse_byte_mask64(it.chunk, '\n') & (~UINT64_C(0) << (uintptr_t)(beg - it.chunk)) : 0;
    return it;
}

MD_PARSE_SIMD_NO_SANITIZE static inline bool md_parse_lines_next(md_parse_lines_t* it, md_str_t* line) {
    const char* nl = it->end;
    for (;;) {
        if (it->mask) {
            nl = it->chunk + md_ctz64(it->mask);
            it->mask &= it->mask - 1;
            if (nl > it->end) nl = it->end;
            break;
        }
        it->chunk += 64;
        if (it->chunk >= it->end) break;
        it->mask = md_parse_byte_mask64(it->chunk, '\n');
    }
    if (it->beg >= it->end) return false;
    const char* line_end = (nl > it->beg && nl[-1] == '\r') ? nl - 1 : nl;
    *line = (md_str_t){ it->beg, (size_t)(line_end - it->beg) };
    it->beg = nl < it->end ? nl + 1 : it->end;
    if (it->beg >= it->end) it->mask = 0;
    return true;
}

// Two fields of 8 characters with 'decimals' decimals (3 or 4), one in each 64 bit lane of v, byte i its i:th
// character: not zero in the lane of one laid out otherwise than md_parse_fixed_f32 takes it (blanks, then an
// optional '-', then digits up to the '.', then the decimals). Of each, its two halves (*parts): the digits of its
// first and second 4 characters as numbers (the '.' a 0), and all bits of its lane where it has no '-' (*pos).
static MD_FORCE_INLINE md_128i md_parse_fixed_w8_pair(md_128i v, int decimals, md_128i* parts, md_128i* pos) {
    const int dot = 7 - decimals;   // 4 or 3
    const md_128i d     = md_mm_sub_epi8(v, md_mm_set1_epi8('0'));
    const md_128i digit = md_mm_cmpeq_epi8(md_mm_min_epu8(d, md_mm_set1_epi8(9)), d);
    const md_128i blank = md_mm_cmpeq_epi8(v, md_mm_set1_epi8(' '));
    const md_128i minus = md_mm_cmpeq_epi8(v, md_mm_set1_epi8('-'));
    const md_128i point = md_mm_cmpeq_epi8(v, md_mm_set1_epi8('.'));
    // The decimals and the last digit before the dot, and the dot
    const uint64_t below_dot = (UINT64_C(1) << (8 * dot)) - 1;
    const uint64_t at_dot    = UINT64_C(0xFF) << (8 * dot);
    const uint64_t want_d    = ~(below_dot | at_dot) | (UINT64_C(0xFF) << (8 * (dot - 1)));
    const md_128i want_digit = md_mm_set1_epi64((int64_t)want_d);
    const md_128i want_dot   = md_mm_set1_epi64((int64_t)at_dot);
    md_128i bad = md_mm_xor_si128(md_mm_and_si128(digit, want_digit), want_digit);
    bad = md_mm_or_si128(bad, md_mm_xor_si128(md_mm_and_si128(point, want_dot), want_dot));
    // What is in front of the digits (the bytes before the last digit before the dot): bytes from the first, which are
    // blanks but for a '-' last
    const md_128i lead = md_mm_andnot_si128(digit, md_mm_set1_epi64((int64_t)below_dot));
    const md_128i last = md_mm_xor_si128(lead, md_mm_srli_epi64(lead, 8));
    bad = md_mm_or_si128(bad, md_mm_and_si128(lead, md_mm_add_epi64(lead, md_mm_set1_epi64(1))));
    bad = md_mm_or_si128(bad, md_mm_andnot_si128(md_mm_or_si128(blank, minus), lead));
    bad = md_mm_or_si128(bad, md_mm_andnot_si128(last, md_mm_and_si128(minus, lead)));
    // Pairs of digits (10, 1), then pairs of those (100, 1): of each field the digits of its two halves as numbers. The
    // dot at 4: the integer part and the decimals; at 3: the integer part times 10 (its ones where the dot is) and the
    // decimals. Either way the field's integer is the first times 1000 and the second.
    const md_128i pairs = md_mm_maddubs_epi16(md_mm_and_si128(d, digit), md_mm_set1_epi16(0x010A));
    *parts = md_mm_madd_epi16(pairs, md_mm_set1_epi32(0x00010064));
    *pos = md_mm_cmpeq_epi64(minus, md_mm_setzero_si128());
    return bad;
}

// The three fields of 8 characters at p ("%8.3f" for decimals 3, "%8.4f" for 4), of which all 24 bytes are read: false,
// out untouched, where any is laid out otherwise than md_parse_fixed_f32 takes it. The value is md_parse_fixed_f32's:
// the field's integer (at most 9999999, exact as a float) over 10^decimals, correctly rounded as a float, as is its
// double rounded to a float (no decimal of 7 digits is near enough to a float's rounding boundary to round otherwise),
// and negated last (a "-0.000" is -0).
static MD_FORCE_INLINE bool md_parse_fixed_float3_w8(float out[3], const char* p, int decimals) {
    MD_ASSERT(decimals == 3 || decimals == 4);
    const md_128i z = md_mm_loadl_epi64((const md_128i*)(p + 16));
    md_128i parts01, parts2, pos01, pos2;
    const md_128i bad = md_mm_or_si128(md_parse_fixed_w8_pair(md_mm_loadu_si128((const md_128i*)p), decimals, &parts01, &pos01),
                                       md_parse_fixed_w8_pair(md_mm_unpacklo_epi64(z, z), decimals, &parts2, &pos2));     // The third twice
    if (!md_mm_testz_si128(bad, bad)) return false;
    // The two halves (at most 9999) as 16 bits, then each field's integer, (1000, 1): of the three, and the third again
    const md_128i n   = md_mm_madd_epi16(md_mm_packs_epi32(parts01, parts2), md_mm_set1_epi32(0x000103E8));
    // Of the low 32 bits of each 64 bit lane
    const md_128i neg = md_mm_castps_si128(md_mm_shuffle_ps(md_mm_castsi128_ps(pos01), md_mm_castsi128_ps(pos2), MD_SHUFFLE(2, 0, 2, 0)));
    const md_128  val = md_mm_div_ps(md_mm_cvtepi32_ps(n), md_mm_set1_ps(decimals == 3 ? 1000.0f : 10000.0f));
    const md_128i sign = md_mm_andnot_si128(neg, md_mm_set1_epi32((int)0x80000000u));
    float res[4];
    md_mm_storeu_ps(res, md_mm_xor_ps(val, md_mm_castsi128_ps(sign)));
    out[0] = res[0];
    out[1] = res[1];
    out[2] = res[2];
    return true;
}

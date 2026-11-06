#pragma once

// FIXED WIDTH FIELDS, THREE AT ONCE
//
// md_parse_fixed_f32 (md_parse.h) reads one field of a fixed width column format. The coordinates of the structure
// formats come three at a time, in fields of 8 characters: "%8.3f" for those of pdb and gro, "%8.4f" for the velocities
// of gro. md_parse_fixed_f32x3_w8 reads three such fields at once (md_simd: SSE/AVX on x86, NEON on arm64), and takes
// exactly what md_parse_fixed_f32 takes, with the same value: where any of the three is laid out otherwise it gives
// false with out untouched, and the caller reads them one by one.
//
// Its own header, so that the includers of md_parse.h do not take md_simd.h with it.

#include <core/md_common.h>
#include <core/md_simd.h>

#include <stdbool.h>
#include <stdint.h>

// Two fields of 8 characters with 'decimals' decimals (3 or 4), one in each 64 bit lane of v, byte i its i:th
// character: not zero in the lane of one laid out otherwise than md_parse_fixed_f32 takes it (blanks, then an
// optional '-', then digits up to the '.', then the decimals). Of each, its two halves (*parts): the digits of its
// first and second 4 characters as numbers (the '.' a 0), and all bits of its lane where it has no '-' (*pos).
static FORCE_INLINE md_128i md_parse_fixed_w8_pair(md_128i v, int decimals, md_128i* parts, md_128i* pos) {
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
static FORCE_INLINE bool md_parse_fixed_f32x3_w8(float out[3], const char* p, int decimals) {
    ASSERT(decimals == 3 || decimals == 4);
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

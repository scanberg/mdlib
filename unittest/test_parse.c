#include "utest.h"

#include <core/md_parse.h>
#include <core/md_reader.h>
#include <core/md_str.h>
#include <core/md_os.h>
#include <core/md_file.h>
#include <core/md_allocator.h>

#include "comma_locale.h"
#include "test_executor.h"

#include <float.h>
#include <math.h>
#include <locale.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// The library and these tests are built with fast math, under which NAN does not compare as itself and
// infinities may be assumed away. Everything about a double's value is therefore checked on its bits.

static uint64_t bits_of(double d) {
    uint64_t u;
    memcpy(&u, &d, sizeof(u));
    return u;
}

static bool is_nan_bits(double d) {
    return (bits_of(d) & 0x7FFFFFFFFFFFFFFFull) > 0x7FF0000000000000ull;
}

static md_str_t str_of(const char* cstr) {
    md_str_t s = { cstr, strlen(cstr) };
    return s;
}

// Deterministic, so a failure reproduces
static uint64_t rng_next(uint64_t* state) {
    uint64_t x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    return *state = x;
}

// ---------------------------------------------------------------------------
// md_parse_f64
// ---------------------------------------------------------------------------

UTEST(parse, f64_grammar) {
    static const struct {
        const char* in;
        size_t      consumed;
        double      value;
    } cases[] = {
        { "0",          1,  0.0 },
        { "1.5",        3,  1.5 },
        { "+1.5",       4,  1.5 },
        { "-1.5",       4, -1.5 },
        { ".5",         2,  0.5 },
        { "-.5e1",      5, -5.0 },
        { "5.",         2,  5.0 },
        { "00012",      5, 12.0 },
        { "0.000",      5,  0.0 },
        { "1e5",        3,  1e5 },
        { "1E+05",      5,  1e5 },
        { "1e-5",       4,  1e-5 },
        { "2.5e-3x",    6,  2.5e-3 },
        // An exponent marker without digits is not part of the number, as for strtod
        { "1e",         1,  1.0 },
        { "1e+",        1,  1.0 },
        { "1ex",        1,  1.0 },
        // A prefix is read, whatever follows it
        { "12.5abc",    4, 12.5 },
        { "1.2.3",      3,  1.2 },
        { "1,5",        1,  1.0 },
        { "0x10",       1,  0.0 },
        { "--1",        0,  0.0 },
        // Out of range
        { "1e400",      5,  HUGE_VAL },
        { "-1e400",     6, -HUGE_VAL },
        { "1e-400",     6,  0.0 },
        { "1e99999999999", 13, HUGE_VAL },
        { "0e99999999999", 13, 0.0 },
        // Not a number at all
        { "",           0,  0.0 },
        { "-",          0,  0.0 },
        { "+",          0,  0.0 },
        { ".",          0,  0.0 },
        { "-.",         0,  0.0 },
        { "e5",         0,  0.0 },
        { " 1",         0,  0.0 },   // md_parse_f64 does not skip whitespace
        { "in",         0,  0.0 },
        // Special values, any case
        { "inf",        3,  HUGE_VAL },
        { "INF",        3,  HUGE_VAL },
        { "-Infinity",  9, -HUGE_VAL },
        { "+infinity",  9,  HUGE_VAL },
        { "info",       3,  HUGE_VAL },
    };

    for (size_t i = 0; i < MD_COUNTOF(cases); ++i) {
        const double sentinel = 12345.0;
        double v = sentinel;
        const size_t n = md_parse_f64(&v, str_of(cases[i].in));
        EXPECT_EQ(cases[i].consumed, n);
        if (n != cases[i].consumed) {
            printf("  input '%s'\n", cases[i].in);
        }
        if (n == 0) {
            // Nothing read, nothing written
            EXPECT_EQ(bits_of(sentinel), bits_of(v));
        } else {
            EXPECT_EQ(bits_of(cases[i].value), bits_of(v));
            if (bits_of(cases[i].value) != bits_of(v)) {
                printf("  input '%s': got %.17g, expected %.17g\n", cases[i].in, v, cases[i].value);
            }
        }
    }

    // The sign of zero survives
    double v = 0.0;
    EXPECT_EQ(2u, md_parse_f64(&v, MD_STR_LIT("-0")));
    EXPECT_EQ(bits_of(-0.0), bits_of(v));
    EXPECT_EQ(7u, md_parse_f64(&v, MD_STR_LIT("-1e-400")));
    EXPECT_EQ(bits_of(-0.0), bits_of(v));

    EXPECT_EQ(3u, md_parse_f64(&v, MD_STR_LIT("NaN")));
    EXPECT_TRUE(is_nan_bits(v));
    EXPECT_EQ(4u, md_parse_f64(&v, MD_STR_LIT("-nan(ind)")));   // The payload MSVC prints is not read
    EXPECT_TRUE(is_nan_bits(v));

    // A null view is no number
    v = 1.0;
    EXPECT_EQ(0u, md_parse_f64(&v, (md_str_t){0}));
    EXPECT_EQ(bits_of(1.0), bits_of(v));
}

// The view is all there is: what lies past its end in memory is not looked at
UTEST(parse, f64_reads_only_its_view) {
    const char text[] = "123.456e7";
    double v = 0.0;

    EXPECT_EQ(3u, md_parse_f64(&v, (md_str_t){ text, 3 }));
    EXPECT_EQ(123.0, v);
    EXPECT_EQ(4u, md_parse_f64(&v, (md_str_t){ text, 4 }));
    EXPECT_EQ(123.0, v);
    EXPECT_EQ(6u, md_parse_f64(&v, (md_str_t){ text, 6 }));
    EXPECT_EQ(bits_of(123.45), bits_of(v));
    EXPECT_EQ(7u, md_parse_f64(&v, (md_str_t){ text, 8 }));    // "123.456e": the 'e' has no digits in view
    EXPECT_EQ(bits_of(123.456), bits_of(v));
    EXPECT_EQ(9u, md_parse_f64(&v, (md_str_t){ text, 9 }));
    EXPECT_EQ(bits_of(123.456e7), bits_of(v));

    const char inf[] = "infinity";
    EXPECT_EQ(3u, md_parse_f64(&v, (md_str_t){ inf, 5 }));     // "infin" reads as "inf"
    EXPECT_EQ(0u, md_parse_f64(&v, (md_str_t){ inf, 2 }));

    EXPECT_TRUE(md_parse_is_float((md_str_t){ text, 3 }));
    EXPECT_TRUE(md_parse_is_float((md_str_t){ text, 4 }));
    EXPECT_FALSE(md_parse_is_float((md_str_t){ text, 8 }));
    EXPECT_TRUE(md_parse_is_int((md_str_t){ text, 3 }));
    EXPECT_FALSE(md_parse_is_int((md_str_t){ text, 4 }));
}

// Inputs where a parser that is merely close gets the last bit wrong: half way cases, the subnormal
// boundary, the edges of the range and mantissas longer than 19 digits.
UTEST(parse, f64_hard_cases) {
    static const char* inputs[] = {
        "9007199254740993",                         // 2^53 + 1: half way, rounds to even
        "9007199254740995",
        "9007199254740993.0000000000000000001",     // Just above half way: rounds up
        "1.00000000000000011102230246251565404236316680908203125",      // Exactly between 1 and its successor
        "1.00000000000000011102230246251565404236316680908203125000001",
        "1.00000000000000011102230246251565404236316680908203124999999",
        "0.1000000000000000055511151231257827021181583404541015625",    // The double nearest 0.1, exactly
        "3.14159265358979323846264338327950288419716939937510582097494459",
        "2.2250738585072011e-308",                  // Below the smallest normal
        "2.2250738585072012e-308",
        "2.2250738585072014e-308",                  // The smallest normal
        "4.9406564584124654e-324",                  // The smallest subnormal
        "2.4703282292062327e-324",                  // Just below half of it: zero
        "2.4703282292062328e-324",                  // Just above: the smallest subnormal
        "1.7976931348623157e308",                   // The largest double
        "1.7976931348623158e308",
        "1.7976931348623159e308",                   // Overflows
        "8.98846567431158e307",
        "1e22", "1e23", "1e-22", "1e-23",           // The edge of the exact powers of ten
        "123456789012345678901234567890",
        "0.000000000000000000000000000000000000001234567890123456789012345",
        "7.2057594037927933e16",
        "4503599627370496.5",                       // Half way at the edge of the exact integers
        "4503599627370497.5",
    };

    for (size_t i = 0; i < MD_COUNTOF(inputs); ++i) {
        const double ref = strtod(inputs[i], NULL);
        double v = 0.0;
        const size_t n = md_parse_f64(&v, str_of(inputs[i]));
        EXPECT_EQ(strlen(inputs[i]), n);
        EXPECT_EQ(bits_of(ref), bits_of(v));
        if (bits_of(ref) != bits_of(v)) {
            printf("  '%s': got %.17g, strtod %.17g\n", inputs[i], v, ref);
        }
    }
}

// Random doubles, printed every way a file is likely to print them, read back to the same bits as the C
// library reads them
UTEST(parse, f64_matches_strtod) {
    uint64_t state = 0x9E3779B97F4A7C15ull;
    size_t mismatches = 0;
    const size_t count = 200000;
    char buf[128];

    for (size_t i = 0; i < count; ++i) {
        // Every finite double is as likely as any other, which spreads the exponents across the range
        uint64_t b = rng_next(&state) & 0x7FFFFFFFFFFFFFFFull;
        if ((b >> 52) == 0x7FF) continue;
        double d;
        memcpy(&d, &b, sizeof(d));
        if (rng_next(&state) & 1) d = -d;

        switch (i % 5) {
        case 0: snprintf(buf, sizeof(buf), "%.17g", d); break;
        case 1: snprintf(buf, sizeof(buf), "%.*g", (int)(rng_next(&state) % 17) + 1, d); break;
        case 2: snprintf(buf, sizeof(buf), "%.*e", (int)(rng_next(&state) % 25), d); break;
        case 3: snprintf(buf, sizeof(buf), "%.*f", (int)(rng_next(&state) % 9), fmod(d, 1e6)); break;
        case 4: snprintf(buf, sizeof(buf), "%.3f", (double)(int64_t)(rng_next(&state) % 2000001) / 1000.0 - 1000.0); break;
        }

        const double ref = strtod(buf, NULL);
        double v = 0.0;
        const size_t n = md_parse_f64(&v, str_of(buf));
        if (n != strlen(buf) || bits_of(v) != bits_of(ref)) {
            if (mismatches++ < 10) {
                printf("  '%s': read %zu of %zu, got %.17g, strtod %.17g\n", buf, n, strlen(buf), v, ref);
            }
        }
    }
    EXPECT_EQ(0u, mismatches);
}

// The seams between the inline part of md_parse_f64 (md_parse.h), md_parse_f64_slow and md_parse_f64_special: each
// input on both sides of a limit one of them checks, read to the same value and length as strtod reads it
UTEST(parse, f64_inline_boundaries) {
    static const char* inputs[] = {
        // 19 digits (the most the scan accumulates) and 20
        "1234567890123456789", "12345678901234567890", "0.1234567890123456789", "1.234567890123456789",
        // The mantissa at 2^53, one past it and far past it
        "9007199254740992", "9007199254740993", "-9007199254740993e-3", "18446744073709551615",
        // Powers of ten at the ends of the exact range, and past them
        "1e-22", "1e-23", "1.5e-22", "123456789e-22", "1e22", "1e23", "1.5E+03", "1.23456789E+02",
        // Exponents of 1 to 4 digits, with and without a sign, and markers without digits
        "1e5", "1e05", "1e+005", "1e-0005", "1e999", "-1e-999", "1e0400", "1e", "1e+", "1E-", "1ex", "1.e5",
        // Zeros of both signs, with exponents the exact path takes and ones it does not
        "0", "-0", "-0.000", "0e0", "-0e-5", "-0e30", "0e9999", "-.0",
        // %E as prmtop and cube files print it
        " 1.23456789E+00", "-1.23456789E-05", "6.02214076E+23", "-4.21000E-12",
    };
    for (size_t i = 0; i < MD_COUNTOF(inputs); ++i) {
        const char* in = inputs[i];
        char* ref_end = NULL;
        const double ref = strtod(in, &ref_end);
        const size_t ref_len = in[0] == ' ' ? 0 : (size_t)(ref_end - in);   // md_parse_f64 does not skip whitespace
        double v = 42.0;
        const size_t n = md_parse_f64(&v, str_of(in));
        EXPECT_EQ_MSG(ref_len, n, in);
        if (ref_len) {
            EXPECT_EQ_MSG(bits_of(ref), bits_of(v), in);
        }
    }
}

// The value of "1.5" does not depend on the locale, which a GTK file dialog (among others) may switch
// to one with a decimal comma
UTEST(parse, f64_ignores_the_locale) {
    static const char* inputs[] = {
        "1.5",
        "-0.25e2",
        "3.14159265358979323846264338327950288",    // Past 19 digits, down the slow path
        "123456789012345678901234567890.5e-10",
    };
    double ref[MD_COUNTOF(inputs)];
    for (size_t i = 0; i < MD_COUNTOF(inputs); ++i) {
        ref[i] = strtod(inputs[i], NULL);
    }

    comma_locale_t loc;
    if (!comma_locale_enter(&loc)) {
        return;
    }

    for (size_t i = 0; i < MD_COUNTOF(inputs); ++i) {
        double v = 0.0;
        EXPECT_EQ(strlen(inputs[i]), md_parse_f64(&v, str_of(inputs[i])));
        EXPECT_EQ(bits_of(ref[i]), bits_of(v));
    }
    comma_locale_leave(&loc);
}

// ---------------------------------------------------------------------------
// md_parse_i64
// ---------------------------------------------------------------------------

UTEST(parse, i64) {
    static const struct {
        const char* in;
        size_t      consumed;
        int64_t     value;
    } cases[] = {
        { "0",                          1,  0 },
        { "-0",                         2,  0 },
        { "+42",                        3,  42 },
        { "-248",                       4, -248 },
        { "000123",                     6,  123 },
        { "12abc",                      2,  12 },
        { "1.9",                        1,  1 },
        { "9223372036854775807",        19, INT64_MAX },
        { "-9223372036854775808",       20, INT64_MIN },
        // Out of range saturates, and every digit is still read
        { "9223372036854775808",        19, INT64_MAX },
        { "-9223372036854775809",       20, INT64_MIN },
        { "18446744073709551615",       20, INT64_MAX },
        { "18446744073709551616",       20, INT64_MAX },
        { "99999999999999999999999",    23, INT64_MAX },
        { "-99999999999999999999999",   24, INT64_MIN },
    };
    for (size_t i = 0; i < MD_COUNTOF(cases); ++i) {
        int64_t v = 7;
        EXPECT_EQ(cases[i].consumed, md_parse_i64(&v, str_of(cases[i].in)));
        EXPECT_EQ(cases[i].value, v);
    }

    static const char* not_ints[] = { "", "-", "+", " 1", "x1", "--1", ".5" };
    for (size_t i = 0; i < MD_COUNTOF(not_ints); ++i) {
        int64_t v = 7;
        EXPECT_EQ(0u, md_parse_i64(&v, str_of(not_ints[i])));
        EXPECT_EQ(7, v);
    }

    const char text[] = "12345";
    int64_t v = 0;
    EXPECT_EQ(2u, md_parse_i64(&v, (md_str_t){ text, 2 }));
    EXPECT_EQ(12, v);
}

UTEST(parse, i64_matches_strtoll) {
    uint64_t state = 0xD1B54A32D192ED03ull;
    char buf[64];
    for (int i = 0; i < 100000; ++i) {
        const int64_t x = (int64_t)(rng_next(&state) >> (rng_next(&state) % 64));
        snprintf(buf, sizeof(buf), "%lld", (long long)((i & 1) ? -x : x));
        int64_t v = 0;
        ASSERT_EQ(strlen(buf), md_parse_i64(&v, str_of(buf)));
        ASSERT_EQ((int64_t)strtoll(buf, NULL, 10), v);
    }
}

// ---------------------------------------------------------------------------
// The conveniences
// ---------------------------------------------------------------------------

UTEST(parse, parse_float_and_parse_int) {
    EXPECT_EQ(-125.0, md_parse_float(MD_STR_LIT("  -12.5e1  ")));
    EXPECT_EQ(1.5, md_parse_float(MD_STR_LIT("\t+1.5")));
    EXPECT_EQ(0.0, md_parse_float(MD_STR_LIT("abc")));
    EXPECT_EQ(0.0, md_parse_float(MD_STR_LIT("")));
    EXPECT_EQ(0.0, md_parse_float((md_str_t){0}));
    EXPECT_EQ(1023.2231128379817, md_parse_float(MD_STR_LIT("1023.22311283798172389718923789172389")));
    EXPECT_EQ(1232326745e10, md_parse_float(MD_STR_LIT("1232326745e10")));
    EXPECT_EQ(1.0e-29, md_parse_float(MD_STR_LIT("1.0e-29")));
    EXPECT_EQ(0.02e+10, md_parse_float(MD_STR_LIT("0.02e+10")));
    EXPECT_EQ(0.273, md_parse_float(MD_STR_LIT("0000000000000.273")));

    EXPECT_EQ(42, md_parse_int(MD_STR_LIT(" 42")));
    EXPECT_EQ(42, md_parse_int(MD_STR_LIT("+42")));
    EXPECT_EQ(-1, md_parse_int(MD_STR_LIT("-1.9")));
    EXPECT_EQ(0, md_parse_int(MD_STR_LIT("abc")));
    EXPECT_EQ(0, md_parse_int(MD_STR_LIT("")));
    EXPECT_EQ(1232326745, md_parse_int(MD_STR_LIT("1232326745")));
}

UTEST(parse, is_float_and_is_int) {
    static const char* floats[] = { "1", "-1", "+1", "1.", ".5", "-.5", "1e5", "1.5E-3", "007", "0.0", "1e+05" };
    static const char* not_floats[] = { "", "-", "+", ".", "e5", "1e", "1e+", "1.5x", " 1", "1 ", "inf", "nan", "Infinity", "1,5", "--1", "1.2.3", "0x10", "1e5.5" };
    for (size_t i = 0; i < MD_COUNTOF(floats); ++i) {
        EXPECT_TRUE(md_parse_is_float(str_of(floats[i])));
    }
    for (size_t i = 0; i < MD_COUNTOF(not_floats); ++i) {
        EXPECT_FALSE(md_parse_is_float(str_of(not_floats[i])));
        if (md_parse_is_float(str_of(not_floats[i]))) printf("  md_parse_is_float('%s')\n", not_floats[i]);
    }
    EXPECT_FALSE(md_parse_is_float((md_str_t){0}));

    static const char* ints[] = { "0", "-0", "+42", "123456789012345678901234567890" };
    static const char* not_ints[] = { "", "-", "+", "1.0", "1e5", " 1", "1 ", "1-" };
    for (size_t i = 0; i < MD_COUNTOF(ints); ++i) {
        EXPECT_TRUE(md_parse_is_int(str_of(ints[i])));
    }
    for (size_t i = 0; i < MD_COUNTOF(not_ints); ++i) {
        EXPECT_FALSE(md_parse_is_int(str_of(not_ints[i])));
    }
    EXPECT_FALSE(md_parse_is_int((md_str_t){0}));
}

// ---------------------------------------------------------------------------
// Tokens
// ---------------------------------------------------------------------------

UTEST(parse, next_token) {
    md_str_t str = MD_STR_LIT("  a bb\tccc\r\n  ");
    md_str_t tok;
    EXPECT_TRUE(md_parse_next_token(&tok, &str));
    EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("a")));
    EXPECT_TRUE(md_parse_next_token(&tok, &str));
    EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("bb")));
    EXPECT_TRUE(md_parse_next_token(&tok, &str));
    EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("ccc")));
    EXPECT_FALSE(md_parse_next_token(&tok, &str));

    md_str_t toks[2];
    str = MD_STR_LIT("1 2 3");
    EXPECT_EQ(2u, md_parse_next_tokens(toks, MD_COUNTOF(toks), &str));
    EXPECT_TRUE(md_str_eq(toks[1], MD_STR_LIT("2")));
    EXPECT_TRUE(md_str_eq(str, MD_STR_LIT("3")));
}

UTEST(parse, next_token_delim) {
    md_str_t toks[8];
    md_str_t str = MD_STR_LIT("a,,b");
    EXPECT_EQ(3u, md_parse_next_tokens_delim(toks, MD_COUNTOF(toks), &str, ','));
    EXPECT_TRUE(md_str_eq(toks[0], MD_STR_LIT("a")));
    EXPECT_TRUE(md_str_eq(toks[1], MD_STR_LIT("")));
    EXPECT_TRUE(md_str_eq(toks[2], MD_STR_LIT("b")));

    // Whitespace is part of a field
    str = MD_STR_LIT(" a , b");
    EXPECT_EQ(2u, md_parse_next_tokens_delim(toks, MD_COUNTOF(toks), &str, ','));
    EXPECT_TRUE(md_str_eq(toks[0], MD_STR_LIT(" a ")));
    EXPECT_TRUE(md_str_eq(toks[1], MD_STR_LIT(" b")));

    // A leading delimiter gives an empty first field; a trailing one no empty last field
    str = MD_STR_LIT(",a,");
    EXPECT_EQ(2u, md_parse_next_tokens_delim(toks, MD_COUNTOF(toks), &str, ','));
    EXPECT_TRUE(md_str_eq(toks[0], MD_STR_LIT("")));
    EXPECT_TRUE(md_str_eq(toks[1], MD_STR_LIT("a")));

    str = MD_STR_LIT("");
    EXPECT_EQ(0u, md_parse_next_tokens_delim(toks, MD_COUNTOF(toks), &str, ','));
}

// ---------------------------------------------------------------------------
// Lines
// ---------------------------------------------------------------------------

UTEST(parse, peek_and_skip_agree_with_extract) {
    md_str_t str = MD_STR_LIT("one\r\ntwo\nlast");
    md_str_t line;

    EXPECT_TRUE(md_parse_peek_line(&line, &str));
    EXPECT_TRUE(md_str_eq(line, MD_STR_LIT("one")));
    EXPECT_TRUE(md_parse_skip_line(&str));
    EXPECT_TRUE(md_parse_peek_line(&line, &str));
    EXPECT_TRUE(md_str_eq(line, MD_STR_LIT("two")));
    EXPECT_TRUE(md_parse_next_line(&line, &str));
    EXPECT_TRUE(md_str_eq(line, MD_STR_LIT("two")));

    // A last line without a newline is still a line, to all three
    EXPECT_TRUE(md_parse_peek_line(&line, &str));
    EXPECT_TRUE(md_str_eq(line, MD_STR_LIT("last")));
    EXPECT_TRUE(md_parse_skip_line(&str));
    EXPECT_FALSE(md_parse_peek_line(&line, &str));
    EXPECT_FALSE(md_parse_skip_line(&str));
    EXPECT_FALSE(md_parse_next_line(&line, &str));
}

// The same lines, through a file buffer that holds only a few of them at a time and through a string
UTEST(parse, reader) {
    md_allocator_i* heap = md_get_heap_allocator();
    const md_str_t path = MD_STR_LIT("md_unittest_parse_lines.txt");

    // Lines of varied length, some CRLF terminated, the last one unterminated
    char* text = md_alloc(heap, 64 * 1024);
    size_t len = 0;
    uint64_t state = 42;
    const int num_lines = 500;
    for (int i = 0; i < num_lines; ++i) {
        const int n = (int)(rng_next(&state) % 40);
        len += snprintf(text + len, 64 * 1024 - len, "%d:", i);
        for (int k = 0; k < n; ++k) text[len++] = (char)('a' + k % 26);
        if (i + 1 < num_lines) {
            if (i % 7 == 0) text[len++] = '\r';
            text[len++] = '\n';
        }
    }

    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
    ASSERT_EQ(len, md_file_write(file, text, len));
    md_file_close(&file);

    // A buffer which fits a few lines at a time, and never all of them
    md_file_reader_t frd;
    ASSERT_TRUE(md_file_reader_open(&frd, path, 64, heap));
    md_reader_i* fr = &frd.r;
    md_reader_i sr = md_reader_from_str((md_str_t){ text, len });
    EXPECT_EQ((int64_t)len, md_reader_size(fr));
    EXPECT_EQ((int64_t)len, md_reader_size(&sr));

    int64_t offset[500];
    for (int pass = 0; pass < 2; ++pass) {
        int count = 0;
        md_str_t a, b, peek;
        while (true) {
            const int64_t pos_f = md_reader_tell(fr);
            const int64_t pos_s = md_reader_tell(&sr);
            EXPECT_EQ(pos_s, pos_f);

            const bool peeked = md_reader_peek_line(fr, &peek);
            const bool got_f = md_reader_extract_line(fr, &a);
            const bool got_s = md_reader_extract_line(&sr, &b);
            EXPECT_EQ(got_s, got_f);
            EXPECT_EQ(got_f, peeked);
            if (!got_f || !got_s) break;
            EXPECT_TRUE(md_str_eq(a, b));
            EXPECT_TRUE(md_str_eq(a, peek));

            int64_t idx = -1;
            EXPECT_TRUE(md_parse_i64(&idx, a) > 0);
            EXPECT_EQ(count, idx);
            if (count < num_lines) offset[count] = pos_f;
            count += 1;
        }
        EXPECT_EQ(num_lines, count);
        EXPECT_TRUE(md_reader_end(fr));
        EXPECT_TRUE(md_reader_end(&sr));

        EXPECT_TRUE(md_reader_rewind(fr));
        EXPECT_TRUE(md_reader_rewind(&sr));
    }

    // Skipping goes past a line without handing it out
    EXPECT_TRUE(md_reader_skip_line(fr));
    md_str_t line;
    EXPECT_TRUE(md_reader_extract_line(fr, &line));
    EXPECT_EQ(0, strncmp(line.ptr, "1:", 2));

    // An offset taken before a line is the line's: back to it, in either direction, outside what is buffered
    for (int i = num_lines - 1; i >= 0; i -= 37) {
        EXPECT_TRUE(md_reader_seek(fr, offset[i]));
        EXPECT_TRUE(md_reader_seek(&sr, offset[i]));
        md_str_t a, b;
        EXPECT_TRUE(md_reader_extract_line(fr, &a));
        EXPECT_TRUE(md_reader_extract_line(&sr, &b));
        EXPECT_TRUE(md_str_eq(a, b));
        int64_t idx = -1;
        EXPECT_TRUE(md_parse_i64(&idx, a) > 0);
        EXPECT_EQ(i, idx);
    }
    EXPECT_FALSE(md_reader_seek(fr, (int64_t)len + 1));

    // The blocks are whole lines, which put together are the text
    EXPECT_TRUE(md_reader_rewind(fr));
    size_t at = 0;
    bool same = true;
    for (md_str_t block; (block = md_reader_block(fr)).len > 0; ) {
        if (at + block.len < len) {     // all but the last, which is the text's unterminated last line
            same &= block.ptr[block.len - 1] == '\n';
        }
        same &= at + block.len <= len && memcmp(text + at, block.ptr, block.len) == 0;
        at += block.len;
        md_reader_advance(fr, block.len);
    }
    EXPECT_TRUE(same);
    EXPECT_EQ(len, at);

    // And read over threads, at offsets: the same text, also where a block is taken in part (what is left is read again)
    {
        test_executor_t threads = test_executor_threads(3);
        EXPECT_TRUE(md_reader_rewind(fr));
        at = 0;
        same = true;
        size_t num_blocks = 0;
        for (md_str_t block; (block = md_reader_block_parallel(fr, test_executor(&threads))).len > 0; ++num_blocks) {
            // Every other block in part: up to a line in its middle
            size_t take = block.len;
            if (num_blocks % 2 == 0 && block.len > 2) {
                const char* nl = memchr(block.ptr + block.len / 2, '\n', block.len - block.len / 2);
                if (nl) take = (size_t)(nl - block.ptr) + 1;
            }
            same &= at + take <= len && memcmp(text + at, block.ptr, take) == 0;
            at += take;
            md_reader_advance(fr, take);
        }
        EXPECT_TRUE(same);
        EXPECT_EQ(len, at);
        EXPECT_TRUE(md_reader_end(fr));
        EXPECT_FALSE(md_reader_failed(fr));
        // Where the next sequential read begins is where the parallel reads ended
        EXPECT_TRUE(md_reader_rewind(fr));
        md_str_t block = md_reader_block_parallel(fr, test_executor(&threads));
        md_reader_advance(fr, block.len);
        md_str_t line;
        if (md_reader_extract_line(fr, &line)) {
            EXPECT_EQ(0, memcmp(text + block.len, line.ptr, line.len));
        }
    }
    md_file_reader_close(&frd);

    // A line longer than the buffer comes in pieces of it, which put together are the line
    {
        ASSERT_TRUE(md_file_reader_open(&frd, path, 8, heap));
        size_t got = 0;
        bool match = true;
        for (md_str_t piece; md_reader_extract_line(&frd.r, &piece); ) {
            match &= piece.len <= 8;
            // what comes back is the text without its line ends: find it in order
            while (got < len && (text[got] == '\n' || text[got] == '\r')) ++got;
            match &= got + piece.len <= len && memcmp(text + got, piece.ptr, piece.len) == 0;
            got += piece.len;
        }
        EXPECT_TRUE(match);
        md_file_reader_close(&frd);
    }

    remove(path.ptr);
    md_free(heap, text, 64 * 1024);
}

// ---------------------------------------------------------------------------
// Timings, for reference: what the C library costs for the same work
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// md_parse_fixed_f32
// ---------------------------------------------------------------------------

static uint32_t bits_of_f32(float f) {
    uint32_t u;
    memcpy(&u, &f, sizeof(u));
    return u;
}

// What the fixed width path must agree with: the float of the correctly rounded double
static float general_f32(md_str_t field) {
    return (float)md_parse_float(field);
}

// Formats m / 10^decimals right aligned in 'width' characters, as "%*.*f" does, without printf
// (whose decimal point follows the locale). False when it does not fit.
static bool format_fixed(char* out, size_t width, int64_t m, size_t decimals) {
    char tmp[32];
    size_t n = 0;
    const bool negative = m < 0;
    uint64_t u = negative ? (uint64_t)(-m) : (uint64_t)m;
    for (size_t i = 0; i < decimals; ++i) {
        tmp[n++] = (char)('0' + u % 10);
        u /= 10;
    }
    tmp[n++] = '.';
    do {
        tmp[n++] = (char)('0' + u % 10);
        u /= 10;
    } while (u);
    if (negative) tmp[n++] = '-';
    if (n > width) return false;
    for (size_t i = 0; i < width - n; ++i) out[i] = ' ';
    for (size_t i = 0; i < n; ++i) out[width - n + i] = tmp[n - 1 - i];
    return true;
}

// Every value "%8.3f" can hold, the layout of gro and pdb coordinates
UTEST(parse, fixed_f32_every_8_3_field) {
    char field[8];
    size_t mismatches = 0, checked = 0;
    for (int64_t m = -999999; m <= 9999999; ++m) {
        ASSERT_TRUE(format_fixed(field, 8, m, 3));
        const md_str_t s = { field, 8 };
        float v = -1.0f;
        if (!md_parse_fixed_f32(&v, s, 3) || (bits_of_f32(v) != bits_of_f32(general_f32(s)) && m != 0)) {
            if (mismatches++ < 5) printf("  '%.8s'\n", field);
        }
        checked += 1;
    }
    EXPECT_EQ(0, mismatches);
    EXPECT_EQ(10999999, checked);
}

// Every width up to 16 and every number of decimals it has room for, on random values, read from an
// exact size copy so a read outside the field is a read outside the allocation
UTEST(parse, fixed_f32_layouts) {
    uint64_t rng = 0x5DEECE66Dull;
    size_t mismatches = 0;
    for (size_t width = 3; width <= 16; ++width) {
        for (size_t decimals = 1; decimals + 2 <= width; ++decimals) {
            for (int k = 0; k < 2000; ++k) {
                // Any number of digits the width holds, so short and long values both
                const size_t digits = 1 + rng_next(&rng) % (width - 1);
                uint64_t limit = 1;
                for (size_t d = 0; d < digits && d < 18; ++d) limit *= 10;
                int64_t m = (int64_t)(rng_next(&rng) % limit);
                if (rng_next(&rng) & 1) m = -m;
                char tmp[32];
                if (!format_fixed(tmp, width, m, decimals)) continue;

                char* field = malloc(width);
                memcpy(field, tmp, width);
                const md_str_t s = { field, width };
                float v = -1.0f;
                const bool ok = md_parse_fixed_f32(&v, s, decimals);
                if (!ok || (bits_of_f32(v) != bits_of_f32(general_f32(s)) && m != 0)) {
                    if (mismatches++ < 5) printf("  '%.*s' with %zu decimals\n", (int)width, tmp, decimals);
                }
                free(field);
            }
        }
    }
    EXPECT_EQ(0, mismatches);
}

// Any other layout is false, with nothing written, however much of a number it is
UTEST(parse, fixed_f32_other_layouts) {
    static const struct {
        const char* field;
        size_t      decimals;
        bool        ok;
        float       value;
    } cases[] = {
        { "   1.234", 3, true,  1.234f },
        { "  -1.234", 3, true, -1.234f },
        { "-123.456", 3, true, -123.456f },
        { "1234.567", 3, true,  1234.567f },
        { "   0.000", 3, true,  0.0f },
        { "  1.00",   2, true,  1.0f },
        { "0.5",      1, true,  0.5f },
        { "-1234.5678", 4, true, -1234.5678f },
        { "  12345.678901", 6, true, 12345.678901f },
        // Not the decimals asked for
        { "  1.2345", 3, false, 0 },
        { "   1.23 ", 3, false, 0 },
        { " 1.23456", 3, false, 0 },
        // No digit in front of the point, or none at all
        { "    .234", 3, false, 0 },
        { "   -.234", 3, false, 0 },
        { "        ", 3, false, 0 },
        { "     ...", 3, false, 0 },
        // Signs: one '-', last, and no '+'
        { "  +1.234", 3, false, 0 },
        { " --1.234", 3, false, 0 },
        { " - 1.234", 3, false, 0 },
        { " -1-.234", 3, false, 0 },
        { "  1-.234", 3, false, 0 },
        // Anything but spaces in front
        { " 1 2.345", 3, false, 0 },
        { "\t  1.234", 3, false, 0 },
        { "  01.234", 3, true,  1.234f },     // leading zeros are digits
        { "x  1.234", 3, false, 0 },
        { "   1.2a4", 3, false, 0 },
        { "   1.23 ", 2, false, 0 },
        { "  1e+03 ", 3, false, 0 },
        { "\xC3\xA9 1.234", 3, false, 0 },
        { "  1.234\xC3", 3, false, 0 },
        { "12345678", 3, false, 0 },
        // Wider than a word, the same rules
        { "    1.23456", 5, true,  1.23456f },
        { "    1.2345 ", 5, false, 0 },
        { "  - 1.23456", 5, false, 0 },
        { "   -.123456", 6, false, 0 },
        { "  1 1.23456", 5, false, 0 },
        // Nothing to read
        { "",         3, false, 0 },
        { "1.5",      0, false, 0 },
        { ".5",       1, false, 0 },
    };
    for (size_t i = 0; i < MD_COUNTOF(cases); ++i) {
        const size_t len = strlen(cases[i].field);
        char* field = malloc(len ? len : 1);
        memcpy(field, cases[i].field, len);
        float v = 42.0f;
        const bool ok = md_parse_fixed_f32(&v, (md_str_t){ field, len }, cases[i].decimals);
        EXPECT_EQ_MSG(cases[i].ok, ok, cases[i].field);
        if (cases[i].ok) {
            EXPECT_EQ_MSG(bits_of_f32(general_f32((md_str_t){ field, len })), bits_of_f32(v), cases[i].field);
            EXPECT_EQ_MSG(cases[i].value, v, cases[i].field);
        } else {
            EXPECT_EQ_MSG(bits_of_f32(42.0f), bits_of_f32(v), cases[i].field);
        }
        free(field);
    }
}

UTEST(parse, perf) {
    static const char* floats[] = { "-248.271233", "12.345", "4.370019985270747", "1.396125e-30", "1232326745" };
    static const char* ints[]   = { "128326746123", "-42", "7" };
    const int num_iter = 200000;
    double facc = 0;
    int64_t iacc = 0;

    md_tick_t t0 = md_tick_now();
    for (int i = 0; i < num_iter; ++i) for (size_t k = 0; k < MD_COUNTOF(floats); ++k) facc += strtod(floats[k], NULL);
    md_tick_t t1 = md_tick_now();
    for (int i = 0; i < num_iter; ++i) for (size_t k = 0; k < MD_COUNTOF(floats); ++k) { double v; md_parse_f64(&v, str_of(floats[k])); facc += v; }
    md_tick_t t2 = md_tick_now();
    for (int i = 0; i < num_iter; ++i) for (size_t k = 0; k < MD_COUNTOF(ints); ++k) iacc += strtoll(ints[k], NULL, 10);
    md_tick_t t3 = md_tick_now();
    for (int i = 0; i < num_iter; ++i) for (size_t k = 0; k < MD_COUNTOF(ints); ++k) { int64_t v; md_parse_i64(&v, str_of(ints[k])); iacc += v; }
    md_tick_t t4 = md_tick_now();

    const double nf = (double)num_iter * MD_COUNTOF(floats);
    const double ni = (double)num_iter * MD_COUNTOF(ints);
    printf("  per number: strtod %.1f ns, md_parse_f64 %.1f ns | strtoll %.1f ns, md_parse_i64 %.1f ns (%g %lld)\n",
        md_tick_to_milliseconds(t1 - t0) * 1e6 / nf, md_tick_to_milliseconds(t2 - t1) * 1e6 / nf,
        md_tick_to_milliseconds(t3 - t2) * 1e6 / ni, md_tick_to_milliseconds(t4 - t3) * 1e6 / ni, facc, (long long)iacc);
}

// The lines of a range in bulk (md_text_lines_t) are those of a plain split, for every start and end in a text of
// empty lines, CRLF ends, long lines, and none at its end: the masks of the chunks before the start and past the end
// see bytes of the text which are not of the range.
static size_t text_lines_reference(md_str_t* out, size_t cap, const char* c, const char* e) {
    size_t n = 0;
    while (c < e && n < cap) {
        const char* nl = c;
        while (nl < e && *nl != '\n') ++nl;
        const char* le = nl;
        if (le > c && le[-1] == '\r') --le;
        out[n++] = (md_str_t){ c, (size_t)(le - c) };
        c = nl < e ? nl + 1 : e;
    }
    return n;
}

UTEST(parse, text_lines) {
    char text[512 + 64];
    uint64_t state = 7;
    for (size_t i = 0; i < 512; ++i) {
        const uint64_t r = rng_next(&state) % 23;
        text[i] = r == 0 ? '\n' : r == 1 ? '\r' : (char)('a' + r);
    }
    // A long line, and runs of empty lines
    for (size_t i = 100; i < 230; ++i) text[i] = 'x';
    for (size_t i = 300; i < 310; ++i) text[i] = '\n';
    md_str_t want[600], got[600];
    size_t mismatches = 0;
    for (size_t beg = 0; beg < 200; beg += 1) {
        for (size_t end = beg; end <= 512; end += 7) {
            const size_t n = text_lines_reference(want, MD_COUNTOF(want), text + beg, text + end);
            md_text_lines_t it = md_text_lines(text + beg, text + end);
            size_t m = 0;
            md_str_t line;
            while (m < MD_COUNTOF(got) && md_text_lines_next(&it, &line)) got[m++] = line;
            if (m != n) { mismatches += 1; continue; }
            for (size_t k = 0; k < n; ++k) {
                if (got[k].ptr != want[k].ptr || got[k].len != want[k].len) { mismatches += 1; break; }
            }
        }
    }
    EXPECT_EQ((size_t)0, mismatches);
}

// Skipping lines lands where extracting them would: by the chunks of what is read, through a buffer that holds a few
// lines at a time and through a string, of any count, the last line without its newline, past the end
UTEST(parse, reader_skip_lines) {
    md_allocator_i* heap = md_get_heap_allocator();
    const md_str_t path = MD_STR_LIT("md_unittest_skip_lines.txt");

    char* text = md_alloc(heap, 256 * 1024);
    size_t len = 0;
    uint64_t state = 7;
    const int num_lines = 3000;
    for (int i = 0; i < num_lines; ++i) {
        const int n = (int)(rng_next(&state) % 150);
        len += snprintf(text + len, 256 * 1024 - len, "%d:", i);
        for (int k = 0; k < n; ++k) text[len++] = (char)('a' + k % 26);
        if (i + 1 < num_lines) {
            if (i % 5 == 0) text[len++] = '\r';
            text[len++] = '\n';
        }
    }

    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
    ASSERT_EQ(len, md_file_write(file, text, len));
    md_file_close(&file);

    // The newlines of the text directly, from every start, against memchr
    {
        bool same = true;
        for (size_t beg = 0; beg < 300; ++beg) {
            for (size_t n = 0; n < 40; n += 3) {
                size_t left = n;
                const char* got = md_text_skip_newlines(text + beg, text + len, &left);
                const char* ref = text + beg;
                size_t ref_left = n;
                while (ref_left) {
                    const char* nl = memchr(ref, '\n', (size_t)(text + len - ref));
                    if (!nl) { ref = text + len; break; }
                    ref = nl + 1;
                    ref_left -= 1;
                }
                same &= got == ref && left == ref_left;
            }
        }
        EXPECT_TRUE(same);
    }

    const size_t caps[] = { 256, 4096, MD_MEGABYTES(1) };
    for (size_t c = 0; c < MD_COUNTOF(caps); ++c) {
        md_file_reader_t frd;
        ASSERT_TRUE(md_file_reader_open(&frd, path, caps[c], heap));
        md_reader_i sr = md_reader_from_str((md_str_t){ text, len });
        md_reader_i* readers[2] = { &frd.r, &sr };
        for (int r = 0; r < 2; ++r) {
            uint64_t st = 11;
            int at = 0;
            bool same = true;
            while (at < num_lines) {
                const size_t k = (size_t)(rng_next(&st) % 64);
                const size_t skipped = md_reader_skip_lines(readers[r], k);
                same &= skipped == MD_MIN(k, (size_t)(num_lines - at));
                at += (int)skipped;
                md_str_t line;
                if (at < num_lines) {
                    same &= md_reader_extract_line(readers[r], &line);
                    int idx = -1;
                    sscanf(line.ptr, "%d:", &idx);
                    same &= idx == at;
                    at += 1;
                }
            }
            EXPECT_TRUE(same);
            EXPECT_EQ(0u, md_reader_skip_lines(readers[r], 5));
            EXPECT_TRUE(md_reader_end(readers[r]));
        }
        md_file_reader_close(&frd);
    }

    remove(path.ptr);
    md_free(heap, text, 256 * 1024);
}

UTEST(parse, next_line) {
    md_str_t str = MD_STR_INIT(
        "this is some text\n"
        "this is line 2\n"
        "\n"
        "\r\n"
        "}\n"
        "this is the end"
    );

    md_str_t line;
    
    EXPECT_TRUE(md_parse_next_line(&line, &str));
    EXPECT_STRNEQ("this is some text", line.ptr, line.len);
    
    EXPECT_TRUE(md_parse_next_line(&line, &str));
    EXPECT_STRNEQ("this is line 2", line.ptr, line.len);
    
    EXPECT_TRUE(md_parse_next_line(&line, &str));
    EXPECT_STRNEQ("", line.ptr, line.len);
    
    EXPECT_TRUE(md_parse_next_line(&line, &str));
    EXPECT_STRNEQ("", line.ptr, line.len);

    EXPECT_TRUE(md_parse_next_line(&line, &str));
    EXPECT_STRNEQ("}", line.ptr, line.len);

    EXPECT_TRUE(md_parse_next_line(&line, &str));
    EXPECT_STRNEQ("this is the end", line.ptr, line.len);

    EXPECT_FALSE(md_parse_next_line(&line, &str));
}

UTEST(parse, skip_peek_line) {
    md_str_t str = MD_STR_LIT("first\r\nsecond\nthird");
    md_str_t line;

    EXPECT_TRUE(md_parse_peek_line(&line, &str));
    EXPECT_TRUE(md_str_eq_cstr(line, "first"));
    EXPECT_TRUE(md_parse_skip_line(&str));
    EXPECT_TRUE(md_parse_next_line(&line, &str));
    EXPECT_TRUE(md_str_eq_cstr(line, "second"));
    EXPECT_TRUE(md_parse_peek_line(&line, &str));
    EXPECT_TRUE(md_str_eq_cstr(line, "third"));
    EXPECT_TRUE(md_parse_skip_line(&str));
    EXPECT_FALSE(md_parse_skip_line(&str));
    EXPECT_FALSE(md_parse_peek_line(&line, &str));

    // A '\r' is only taken off together with a '\n'
    str = MD_STR_LIT("end\r");
    EXPECT_TRUE(md_parse_next_line(&line, &str));
    EXPECT_TRUE(md_str_eq_cstr(line, "end\r"));
}

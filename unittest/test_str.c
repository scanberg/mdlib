#include "utest.h"

#include <core/md_str.h>
#include <core/md_format.h>
#include <core/md_os.h>
#include <core/md_allocator.h>
#include <core/md_writer.h>
#include <md_types.h>
#include <core/md_parse.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include "comma_locale.h"

// What mdlib formats reads back with md_parse, which does not follow the locale either: under a decimal
// comma locale 1.5 is still written "1.5", by md_snprintf and by everything formatting through it
UTEST(str, printf_ignores_the_locale) {
    comma_locale_t loc;
    if (!comma_locale_enter(&loc)) {
        return;
    }

    char buf[128];
    // The C library follows the locale, or this would check nothing
    snprintf(buf, sizeof(buf), "%.2f", 1.5);
    EXPECT_STREQ("1,50", buf);

    EXPECT_EQ(4, md_snprintf(buf, sizeof(buf), "%.2f", 1.5));
    EXPECT_STREQ("1.50", buf);
    md_snprintf(buf, sizeof(buf), "%.9g|%.17g|%e|%g", (double)0.1f, 0.25, 1.5e-7, -2.5);
    EXPECT_STREQ("0.100000001|0.25|1.500000e-07|-2.5", buf);

    // Measured, and cut short, the way C99 has it
    EXPECT_EQ(4, md_snprintf(NULL, 0, "%.2f", 1.5));
    EXPECT_EQ(4, md_snprintf(buf, 3, "%.2f", 1.5));
    EXPECT_STREQ("1.", buf);

    md_allocator_i* alloc = md_get_heap_allocator();
    md_str_t str = md_alloc_printf(alloc, "%.3f", 2.25);
    EXPECT_TRUE(md_str_eq_cstr(str, "2.250"));
    md_free_str(alloc, str);

    md_str_writer_t sw = md_str_writer_create(alloc);
    md_str_writer_fmt(&sw, "%g,%g", 0.5, 1.5);
    EXPECT_TRUE(md_str_eq_cstr(md_str_writer_to_str(&sw), "0.5,1.5"));
    md_str_writer_free(&sw);

    comma_locale_leave(&loc);
}

UTEST(str, printf_edges) {
    char buf[8] = "xxxxxxx";
    // No room at all: measured, and nothing written - not even a terminator in front of the buffer
    EXPECT_EQ(5, md_snprintf(buf + 1, 0, "%s", "hello"));
    EXPECT_STREQ("xxxxxxx", buf);

    // The C libraries spell them differently (glibc "nan", the UCRT "-nan(ind)"); md_parse reads every one
    const double specials[] = { NAN, -NAN, INFINITY, -INFINITY };
    for (size_t i = 0; i < sizeof(specials) / sizeof(specials[0]); ++i) {
        char out[64];
        const int len = md_snprintf(out, sizeof(out), "%f", specials[i]);
        double v = 0;
        EXPECT_NE(0u, md_parse_f64(&v, (md_str_t){out, (size_t)len}));
        EXPECT_EQ(!!isnan(specials[i]), !!isnan(v));
        EXPECT_EQ(!!isinf(specials[i]), !!isinf(v));
        if (isinf(specials[i])) {
            EXPECT_EQ(specials[i] < 0, v < 0);
        }
    }
}

// What is written as %.17g (double) or %.9g (float) reads back with md_parse as the very same value
UTEST(str, printf_round_trips_through_md_parse) {
    uint64_t s = 0x9E3779B97F4A7C15ull;
    size_t fails_f64 = 0, fails_f32 = 0;
    char buf[64];
    for (int i = 0; i < 200000; ++i) {
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        // Any finite bit pattern, subnormals included
        uint64_t bits = s;
        if (((bits >> 52) & 0x7FF) == 0x7FF) continue;
        double d;
        memcpy(&d, &bits, sizeof(d));

        const int len = md_snprintf(buf, sizeof(buf), "%.17g", d);
        double r = 0;
        if (md_parse_f64(&r, (md_str_t){buf, (size_t)len}) != (size_t)len || memcmp(&r, &d, sizeof(d)) != 0) {
            fails_f64 += 1;
        }

        uint32_t fbits = (uint32_t)(s >> 32);
        if (((fbits >> 23) & 0xFF) == 0xFF) continue;
        float f;
        memcpy(&f, &fbits, sizeof(f));
        const int flen = md_snprintf(buf, sizeof(buf), "%.9g", (double)f);
        double rf = 0;
        md_parse_f64(&rf, (md_str_t){buf, (size_t)flen});
        const float back = (float)rf;
        if (memcmp(&back, &f, sizeof(f)) != 0) {
            fails_f32 += 1;
        }
    }
    EXPECT_EQ(0u, fails_f64);
    EXPECT_EQ(0u, fails_f32);
}

UTEST(str, chars_in_place) {
    char buf[] = "C:\\Data\\Run-1.PDB";
    md_chars_replace(buf, sizeof(buf) - 1, '\\', '/');
    EXPECT_STREQ("C:/Data/Run-1.PDB", buf);
    md_chars_to_lower(buf, sizeof(buf) - 1);
    EXPECT_STREQ("c:/data/run-1.pdb", buf);
    // Only the first len characters, and nothing past them
    md_chars_to_upper(buf, 6);
    EXPECT_STREQ("C:/DATa/run-1.pdb", buf);
    // Only ASCII letters change: UTF-8 bytes are left as they are
    char utf8[] = "\xC3\x85ngstr\xC3\xB6m";
    md_chars_to_upper(utf8, sizeof(utf8) - 1);
    EXPECT_STREQ("\xC3\x85NGSTR\xC3\xB6M", utf8);
    // An empty run, with no storage at all
    md_chars_replace(NULL, 0, 'a', 'b');
    md_chars_to_lower(NULL, 0);
    md_chars_to_upper(NULL, 0);
}

UTEST(str, find_str) {
    EXPECT_EQ((size_t)6, md_str_find_str(MD_STR_LIT("hello world"), MD_STR_LIT("world")));
    EXPECT_EQ((size_t)0, md_str_find_str(MD_STR_LIT("hello"), MD_STR_LIT("hello")));
    EXPECT_EQ((size_t)4, md_str_find_str(MD_STR_LIT("hello"), MD_STR_LIT("o")));

    // The match starts inside an earlier partial match
    EXPECT_EQ((size_t)1, md_str_find_str(MD_STR_LIT("aab"), MD_STR_LIT("ab")));
    EXPECT_EQ((size_t)2, md_str_find_str(MD_STR_LIT("ababac"), MD_STR_LIT("abac")));
    EXPECT_EQ((size_t)3, md_str_find_str(MD_STR_LIT("<a id=\"x\"></a>"), MD_STR_LIT("id=\"")));

    // First occurrence
    EXPECT_EQ((size_t)0, md_str_find_str(MD_STR_LIT("xyxy"), MD_STR_LIT("xy")));

    // Nothing to find
    EXPECT_EQ(MD_STR_NPOS, md_str_find_str(MD_STR_LIT("hello"), MD_STR_LIT("world")));
    EXPECT_EQ(MD_STR_NPOS, md_str_find_str(MD_STR_LIT("hell"),  MD_STR_LIT("hello")));   // needle longer than haystack
    EXPECT_EQ(MD_STR_NPOS, md_str_find_str(MD_STR_LIT("hello"), MD_STR_LIT("")));
    EXPECT_EQ(MD_STR_NPOS, md_str_find_str(MD_STR_LIT(""),      MD_STR_LIT("a")));
    EXPECT_EQ(MD_STR_NPOS, md_str_find_str(MD_STR_LIT("abc"),   MD_STR_LIT("bcd")));     // would run past the end

    EXPECT_TRUE (md_str_contains_str(MD_STR_LIT("abc"), MD_STR_LIT("bc")));
    EXPECT_FALSE(md_str_contains_str(MD_STR_LIT("abc"), MD_STR_LIT("bd")));
}

UTEST(str, copy_to_char_buf) {
    char buf[8];

    EXPECT_EQ(md_str_copy_to_char_buf(buf, sizeof(buf), MD_STR_LIT("abc")), 3);
    EXPECT_STREQ("abc", buf);

    // Truncated to fit, and still zero terminated
    EXPECT_EQ(md_str_copy_to_char_buf(buf, sizeof(buf), MD_STR_LIT("0123456789")), 7);
    EXPECT_STREQ("0123456", buf);

    // An empty string clears the buffer instead of leaving the previous contents
    EXPECT_EQ(md_str_copy_to_char_buf(buf, sizeof(buf), MD_STR_LIT("")), 0);
    EXPECT_STREQ("", buf);
    md_str_copy_to_char_buf(buf, sizeof(buf), MD_STR_LIT("abc"));
    EXPECT_EQ(md_str_copy_to_char_buf(buf, sizeof(buf), (md_str_t){0}), 0);
    EXPECT_STREQ("", buf);

    // A capacity of one only has room for the terminator, a capacity of zero writes nothing
    buf[0] = 'x';
    EXPECT_EQ(md_str_copy_to_char_buf(buf, 1, MD_STR_LIT("abc")), 0);
    EXPECT_EQ(buf[0], '\0');
    buf[0] = 'x';
    EXPECT_EQ(md_str_copy_to_char_buf(buf, 0, MD_STR_LIT("abc")), 0);
    EXPECT_EQ(buf[0], 'x');
}

// A md_str_t is a view: an empty one is the empty string wherever it points, as with std::string_view
UTEST(str, eq_empty) {
    const md_str_t null_str  = {0};
    const md_str_t lit       = MD_STR_LIT("");
    const md_str_t slice     = md_str_substr(MD_STR_LIT("abc"), 1, 0);   // empty, pointing into another string
    const md_str_t slice2    = md_str_substr(MD_STR_LIT("xyz"), 2, 0);

    EXPECT_TRUE(md_str_eq(null_str, null_str));
    EXPECT_TRUE(md_str_eq(null_str, lit));
    EXPECT_TRUE(md_str_eq(lit, null_str));
    EXPECT_TRUE(md_str_eq(slice, slice2));
    EXPECT_TRUE(md_str_eq(slice, null_str));
    EXPECT_FALSE(md_str_eq(null_str, MD_STR_LIT("a")));

    EXPECT_TRUE(md_str_eq_ignore_case(null_str, lit));
    EXPECT_TRUE(md_str_eq_n(null_str, lit, 4));
    EXPECT_TRUE(md_str_eq_n_ignore_case(slice, null_str, 4));

    EXPECT_TRUE(md_str_eq_cstr(null_str, ""));
    EXPECT_TRUE(md_str_eq_cstr(lit, ""));
    EXPECT_TRUE(md_str_eq_cstr(null_str, NULL));          // a NULL cstr is "", as for md_str_from_cstr
    EXPECT_FALSE(md_str_eq_cstr(MD_STR_LIT("a"), NULL));
    EXPECT_FALSE(md_str_eq_cstr(null_str, "a"));
    EXPECT_TRUE(md_str_eq_cstr_ignore_case(slice, ""));
    EXPECT_TRUE(md_str_eq_cstr_n(null_str, "", 3));

    // A copy of an empty string is equal to it ({0} comes back from md_alloc_str)
    md_allocator_i* alloc = md_get_heap_allocator();
    md_str_t copy = md_alloc_str(alloc, lit);
    EXPECT_TRUE(md_str_eq(copy, lit));
    md_free_str(alloc, copy);
}

UTEST(str, eq) {
    EXPECT_TRUE(md_str_eq(MD_STR_LIT("abc"), MD_STR_LIT("abc")));
    EXPECT_FALSE(md_str_eq(MD_STR_LIT("abc"), MD_STR_LIT("abd")));
    EXPECT_FALSE(md_str_eq(MD_STR_LIT("abc"), MD_STR_LIT("ab")));
    EXPECT_TRUE(md_str_eq(md_str_substr(MD_STR_LIT("xabcx"), 1, 3), MD_STR_LIT("abc")));

    EXPECT_TRUE(md_str_eq_ignore_case(MD_STR_LIT("HeLLo"), MD_STR_LIT("hello")));
    EXPECT_FALSE(md_str_eq_ignore_case(MD_STR_LIT("hello"), MD_STR_LIT("hell")));

    EXPECT_TRUE(md_str_eq_cstr(MD_STR_LIT("abc"), "abc"));
    EXPECT_FALSE(md_str_eq_cstr(MD_STR_LIT("abc"), "abcd"));
    EXPECT_FALSE(md_str_eq_cstr(MD_STR_LIT("abc"), "ab"));
    EXPECT_TRUE(md_str_eq_cstr(md_str_substr(MD_STR_LIT("abcdef"), 0, 3), "abc"));   // not zero terminated
    EXPECT_TRUE(md_str_eq_cstr_ignore_case(MD_STR_LIT("ATOM"), "atom"));
    EXPECT_FALSE(md_str_eq_cstr_ignore_case(MD_STR_LIT("ATOM"), "atoms"));
}

// The _n variants compare the first n characters of both: with n = strlen(prefix) they test for a prefix,
// which a string shorter than the prefix does not have
UTEST(str, eq_n) {
    EXPECT_TRUE(md_str_eq_n(MD_STR_LIT("ATOM      1"), MD_STR_LIT("ATOM"), 4));
    EXPECT_TRUE(md_str_eq_n(MD_STR_LIT("abcdef"), MD_STR_LIT("abcxyz"), 3));
    EXPECT_FALSE(md_str_eq_n(MD_STR_LIT("abcdef"), MD_STR_LIT("abcxyz"), 4));
    EXPECT_FALSE(md_str_eq_n(MD_STR_LIT("ab"), MD_STR_LIT("abc"), 3));
    EXPECT_TRUE(md_str_eq_n(MD_STR_LIT("ab"), MD_STR_LIT("ab"), 3));       // both end within n, and are equal
    EXPECT_TRUE(md_str_eq_n(MD_STR_LIT("abc"), MD_STR_LIT("xyz"), 0));

    EXPECT_TRUE(md_str_eq_cstr_n(MD_STR_LIT("ENDMDL"), "ENDMDL", 6));
    EXPECT_TRUE(md_str_eq_cstr_n(MD_STR_LIT("ENDMDL      2"), "ENDMDL", 6));
    EXPECT_FALSE(md_str_eq_cstr_n(MD_STR_LIT("END"), "ENDMDL", 6));
    EXPECT_FALSE(md_str_eq_cstr_n(MD_STR_LIT("C"), "C3", 2));
    EXPECT_TRUE(md_str_eq_cstr_n(MD_STR_LIT("C3'"), "C3", 2));
    EXPECT_TRUE(md_str_eq_cstr_n(MD_STR_LIT("lj/cut"), "lj", 2));
    EXPECT_FALSE(md_str_eq_cstr_n(MD_STR_LIT("REMARK"), "REMARK 350", 10));
    EXPECT_TRUE(md_str_eq_cstr_n(MD_STR_LIT("ab"), "ab", 5));             // both end within n

    EXPECT_TRUE(md_str_eq_n_ignore_case(MD_STR_LIT("Vibration 1"), MD_STR_LIT("VIBRATION"), 9));
    EXPECT_TRUE(md_str_eq_cstr_n_ignore_case(MD_STR_LIT("Time (ps)"), "time", 4));
    EXPECT_FALSE(md_str_eq_cstr_n_ignore_case(MD_STR_LIT("Tim"), "time", 4));
}

UTEST(str, cmp_lex) {
    EXPECT_EQ(md_str_cmp_lex(MD_STR_LIT("abc"), MD_STR_LIT("abc")), 0);
    EXPECT_EQ(md_str_cmp_lex(MD_STR_LIT("abc"), MD_STR_LIT("abd")), -1);
    EXPECT_EQ(md_str_cmp_lex(MD_STR_LIT("abd"), MD_STR_LIT("abc")), 1);
    EXPECT_EQ(md_str_cmp_lex(MD_STR_LIT("ab"),  MD_STR_LIT("abc")), -1);    // a prefix first
    EXPECT_EQ(md_str_cmp_lex(MD_STR_LIT("abc"), MD_STR_LIT("ab")),  1);
    EXPECT_EQ(md_str_cmp_lex((md_str_t){0},     MD_STR_LIT("")),    0);
    EXPECT_EQ(md_str_cmp_lex((md_str_t){0},     MD_STR_LIT("a")),  -1);
    // Characters order as unsigned bytes: UTF-8 (lead byte >= 0x80) after ASCII, whatever the signedness of char
    EXPECT_EQ(md_str_cmp_lex(MD_STR_LIT("z"), MD_STR_LIT("\xC3\xA5")), -1);
}

UTEST(str, begins_ends_with) {
    EXPECT_TRUE(md_str_begins_with(MD_STR_LIT("hello"), MD_STR_LIT("he")));
    EXPECT_FALSE(md_str_begins_with(MD_STR_LIT("he"), MD_STR_LIT("hello")));
    EXPECT_TRUE(md_str_begins_with(MD_STR_LIT("hello"), (md_str_t){0}));
    EXPECT_TRUE(md_str_begins_with((md_str_t){0}, MD_STR_LIT("")));
    EXPECT_TRUE(md_str_ends_with(MD_STR_LIT("file.xtc"), MD_STR_LIT(".xtc")));
    EXPECT_FALSE(md_str_ends_with(MD_STR_LIT("xtc"), MD_STR_LIT(".xtc")));
    EXPECT_TRUE(md_str_ends_with(MD_STR_LIT("file"), (md_str_t){0}));
}

UTEST(str, find_char) {
    EXPECT_EQ((size_t)1, md_str_find_char (MD_STR_LIT("a=b=c"), '='));
    EXPECT_EQ((size_t)3, md_str_rfind_char(MD_STR_LIT("a=b=c"), '='));
    EXPECT_EQ((size_t)0, md_str_rfind_char(MD_STR_LIT("abc"), 'a'));
    EXPECT_EQ(MD_STR_NPOS, md_str_find_char (MD_STR_LIT("abc"), 'x'));
    EXPECT_EQ(MD_STR_NPOS, md_str_rfind_char(MD_STR_LIT("abc"), 'x'));
    EXPECT_EQ(MD_STR_NPOS, md_str_find_char ((md_str_t){0}, 'a'));
    EXPECT_EQ(MD_STR_NPOS, md_str_rfind_char((md_str_t){0}, 'a'));
    EXPECT_TRUE (md_str_contains_char(MD_STR_LIT("abc"), 'c'));
    EXPECT_FALSE(md_str_contains_char(MD_STR_LIT(""), 'c'));

    // A byte above 0x7F is matched as that byte
    const md_str_t utf8 = MD_STR_LIT("x\xC3\xA5y");
    EXPECT_EQ((size_t)1, md_str_find_char (utf8, (char)0xC3));
    EXPECT_EQ((size_t)2, md_str_find_char (utf8, (char)0xA5));
    EXPECT_EQ((size_t)1, md_str_rfind_char(utf8, (char)0xC3));

    // The (ptr, len) layer underneath: only the first len characters, and (NULL, 0) is fine
    EXPECT_EQ((size_t)2,   md_chars_find ("abc", 3, 'c'));
    EXPECT_EQ(MD_STR_NPOS, md_chars_find ("abc", 2, 'c'));
    EXPECT_EQ(MD_STR_NPOS, md_chars_rfind("abc", 2, 'c'));
    EXPECT_EQ(MD_STR_NPOS, md_chars_find (NULL, 0, 'c'));
    EXPECT_EQ(MD_STR_NPOS, md_chars_rfind(NULL, 0, 'c'));
}

// MD_STR_NPOS is SIZE_MAX, which md_str_substr reads as "to the end"
UTEST(str, npos_composes_with_substr) {
    const md_str_t kv = MD_STR_LIT("key=value");
    EXPECT_TRUE(md_str_eq_cstr(md_str_substr(kv, 0, md_str_find_char(kv, '=')), "key"));
    EXPECT_TRUE(md_str_eq_cstr(md_str_substr(kv, 0, md_str_find_char(kv, '#')), "key=value"));

    // After the last '/': MD_STR_NPOS + 1 is 0, which keeps all of it
    const md_str_t path = MD_STR_LIT("a/b/c.txt"), file = MD_STR_LIT("c.txt");
    EXPECT_TRUE(md_str_eq_cstr(md_str_substr(path, md_str_rfind_char(path, '/') + 1, SIZE_MAX), "c.txt"));
    EXPECT_TRUE(md_str_eq_cstr(md_str_substr(file, md_str_rfind_char(file, '/') + 1, SIZE_MAX), "c.txt"));
}

UTEST(str, cut) {
    md_str_t a, b;
    const md_str_t kv = MD_STR_LIT("key=value=x");

    // At the first and at the last '='
    EXPECT_TRUE(md_str_cut(&a, &b, kv, md_str_find_char(kv, '=')));
    EXPECT_TRUE(md_str_eq_cstr(a, "key"));
    EXPECT_TRUE(md_str_eq_cstr(b, "value=x"));
    EXPECT_TRUE(md_str_cut(&a, &b, kv, md_str_rfind_char(kv, '=')));
    EXPECT_TRUE(md_str_eq_cstr(a, "key=value"));
    EXPECT_TRUE(md_str_eq_cstr(b, "x"));

    // At either end, that part is empty
    EXPECT_TRUE(md_str_cut(&a, &b, kv, 0));
    EXPECT_EQ((size_t)0, a.len);
    EXPECT_TRUE(md_str_eq_cstr(b, "ey=value=x"));
    EXPECT_TRUE(md_str_cut(&a, &b, kv, kv.len - 1));
    EXPECT_TRUE(md_str_eq_cstr(a, "key=value="));
    EXPECT_EQ((size_t)0, b.len);

    // Nothing found, or past the end: before is all of it, after is empty
    EXPECT_FALSE(md_str_cut(&a, &b, kv, md_str_find_char(kv, '#')));
    EXPECT_TRUE(md_str_eq_cstr(a, "key=value=x"));
    EXPECT_EQ((size_t)0, b.len);
    EXPECT_FALSE(md_str_cut(&a, &b, kv, kv.len));
    EXPECT_TRUE(md_str_eq_cstr(a, "key=value=x"));
    EXPECT_FALSE(md_str_cut(&a, &b, (md_str_t){0}, 0));
    EXPECT_EQ((size_t)0, a.len);

    // The out pointers may be NULL, and may point at the input: walking a list
    md_str_t list = MD_STR_LIT("a,b,c");
    md_str_t item;
    EXPECT_TRUE(md_str_cut(&item, &list, list, md_str_find_char(list, ',')));
    EXPECT_TRUE(md_str_eq_cstr(item, "a"));
    EXPECT_TRUE(md_str_eq_cstr(list, "b,c"));
    EXPECT_TRUE(md_str_cut(NULL, &list, list, md_str_find_char(list, ',')));
    EXPECT_TRUE(md_str_eq_cstr(list, "c"));
    EXPECT_FALSE(md_str_cut(&item, &list, list, md_str_find_char(list, ',')));
    EXPECT_TRUE(md_str_eq_cstr(item, "c"));
    EXPECT_EQ((size_t)0, list.len);
}

UTEST(str, swap) {
    md_str_t a = MD_STR_LIT("a");
    md_str_t b = MD_STR_LIT("bb");
    md_str_swap(&a, &b);
    EXPECT_TRUE(md_str_eq_cstr(a, "bb"));
    EXPECT_TRUE(md_str_eq_cstr(b, "a"));
}

// A label holds up to 8 characters (md_label16_t 16), zero after them: its length is where the first zero byte is,
// all of it when there is none, and two labels are equal when their bytes are
UTEST(str, labels) {
    static const char text[] = "ABCDEFGHIJKLMNOPQRS";
    for (size_t n = 0; n <= 10; ++n) {
        const md_label_t l = make_label((md_str_t){text, n});
        EXPECT_EQ(MD_MIN(n, 8), md_label_len(&l));
        EXPECT_TRUE(md_str_eq(md_label_str(&l), (md_str_t){text, MD_MIN(n, 8)}));
        EXPECT_EQ(n == 0, label_empty(l));
    }
    for (size_t n = 0; n <= 18; ++n) {
        const md_label16_t l = make_label16((md_str_t){text, n});
        EXPECT_EQ(MD_MIN(n, 16), md_label16_len(&l));
        EXPECT_TRUE(md_str_eq(md_label16_str(&l), (md_str_t){text, MD_MIN(n, 16)}));
    }
    // Bytes past ASCII count as characters
    const md_label_t u = make_label(MD_STR_LIT("C\xc3\xa5\x80\xff"));
    EXPECT_EQ(5, md_label_len(&u));

    EXPECT_TRUE(md_label_eq(make_label(MD_STR_LIT("CA")), make_label(MD_STR_LIT("CA"))));
    EXPECT_FALSE(md_label_eq(make_label(MD_STR_LIT("CA")), make_label(MD_STR_LIT("CA1"))));
    EXPECT_TRUE(md_label_eq(make_label(MD_STR_LIT("ABCDEFGH")), make_label(MD_STR_LIT("ABCDEFGHI"))));
    EXPECT_NE(md_label_key(make_label(MD_STR_LIT("CA"))), md_label_key(make_label(MD_STR_LIT("CB"))));
}

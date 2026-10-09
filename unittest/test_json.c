#include "utest.h"

#include <core/md_json.h>
#include <core/md_parse.h>
#include <core/md_str.h>
#include <core/md_os.h>
#include <core/md_allocator.h>

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Built with fast math like the library, so doubles are compared on their bits

static uint64_t bits_of(double d) {
    uint64_t u;
    memcpy(&u, &d, sizeof(u));
    return u;
}

static str_t str_of(const char* cstr) {
    str_t s = { cstr, strlen(cstr) };
    return s;
}

static uint64_t rng_next(uint64_t* state) {
    uint64_t x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    return *state = x;
}

static uint32_t rng_range(uint64_t* state, uint32_t n) {
    return (uint32_t)(rng_next(state) % n);
}

// Parses from an exact size heap copy of the text, so a read past its end is a read past the
// allocation, which the address sanitizer reports
typedef struct parsed_t {
    char*      text;
    md_json_t* doc;
    md_json_error_t err;
} parsed_t;

static parsed_t parse_n(const char* text, size_t len) {
    parsed_t p = {0};
    p.text = malloc(len ? len : 1);
    if (len) memcpy(p.text, text, len);
    str_t s = { p.text, len };
    p.doc = md_json_parse(s, md_get_heap_allocator(), &p.err);
    return p;
}

static parsed_t parse(const char* text) {
    return parse_n(text, strlen(text));
}

static void release(parsed_t* p) {
    md_json_free(p->doc, md_get_heap_allocator());
    free(p->text);
    p->doc = NULL;
    p->text = NULL;
}

static md_json_val_t get(md_json_val_t v, const char* key) {
    return md_json_get(v, str_of(key));
}

// ---------------------------------------------------------------------------
// Grammar
// ---------------------------------------------------------------------------

UTEST(json, root_values) {
    static const struct {
        const char*    text;
        md_json_type_t type;
    } cases[] = {
        { "null",               MD_JSON_TYPE_NULL   },
        { "true",               MD_JSON_TYPE_BOOL   },
        { "false",              MD_JSON_TYPE_BOOL   },
        { "0",                  MD_JSON_TYPE_NUMBER },
        { "-0",                 MD_JSON_TYPE_NUMBER },
        { "-12.5e-3",           MD_JSON_TYPE_NUMBER },
        { "1E+2",               MD_JSON_TYPE_NUMBER },
        { "NaN",                MD_JSON_TYPE_NUMBER },
        { "Infinity",           MD_JSON_TYPE_NUMBER },
        { "-Infinity",          MD_JSON_TYPE_NUMBER },
        { "\"\"",               MD_JSON_TYPE_STRING },
        { "\"s\"",              MD_JSON_TYPE_STRING },
        { "[]",                 MD_JSON_TYPE_ARRAY  },
        { "{}",                 MD_JSON_TYPE_OBJECT },
        { " \t\r\n 7 \n\t\r ",  MD_JSON_TYPE_NUMBER },
        { "\xEF\xBB\xBF[1]",    MD_JSON_TYPE_ARRAY  },
        { "\xEF\xBB\xBF \"x\"", MD_JSON_TYPE_STRING },
        { "[ ]",                MD_JSON_TYPE_ARRAY  },
        { "{ \n }",             MD_JSON_TYPE_OBJECT },
    };
    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        parsed_t p = parse(cases[i].text);
        EXPECT_TRUE_MSG(p.doc != NULL, cases[i].text);
        EXPECT_EQ_MSG(cases[i].type, md_json_type(md_json_root(p.doc)), cases[i].text);
        // The root has no siblings and no key
        EXPECT_FALSE(md_json_valid(md_json_next(md_json_root(p.doc))));
        EXPECT_FALSE(md_json_valid(md_json_key(md_json_root(p.doc))));
        release(&p);
    }
}

UTEST(json, rejects) {
    static const struct {
        const char* text;
        size_t      offset;     // where the error is reported
    } cases[] = {
        { "",                       0 },
        { "   ",                    3 },
        { "[",                      1 },
        { "[1",                     2 },
        { "[1,",                    3 },
        { "[1,]",                   3 },
        { "[,1]",                   1 },
        { "[1 2]",                  3 },
        { "[1,2,]",                 5 },
        { "]",                      0 },
        { "}",                      0 },
        { "[1]]",                   3 },
        { "[1]x",                   3 },
        { "{}{}",                   2 },
        { "[}",                     1 },
        { "{]",                     1 },
        { "[1}",                    2 },
        { "{\"a\":1]",              6 },
        { "{",                      1 },
        { "{\"a\"",                 4 },
        { "{\"a\":",                5 },
        { "{\"a\":1,",              7 },
        { "{\"a\":1,}",             7 },
        { "{a:1}",                  1 },
        { "{'a':1}",                1 },
        { "{1:1}",                  1 },
        { "{\"a\" 1}",              5 },
        { "{\"a\"::1}",             5 },
        { "{\"a\":1 \"b\":2}",      7 },
        { "{,}",                    1 },
        // Numbers
        { "01",                     0 },
        { "-01",                    0 },
        { "-",                      0 },
        { "-a",                     0 },
        { "+1",                     0 },
        { ".5",                     0 },
        { "1.",                     0 },
        { "1.e5",                   0 },
        { "1e",                     0 },
        { "1e+",                    0 },
        { "1E-",                    0 },
        { "0x10",                   0 },
        { "1.5.2",                  0 },
        { "1a",                     0 },
        { "[1.5x]",                 1 },
        { "- 1",                    0 },
        { "-NaN",                   0 },
        { "-inf",                   0 },
        { "--1",                    0 },
        // Literals
        { "tru",                    0 },
        { "truex",                  0 },
        { "nul",                    0 },
        { "fals",                   0 },
        { "True",                   0 },
        { "NULL",                   0 },
        { "nan",                    0 },
        { "NAN",                    0 },
        { "inf",                    0 },
        { "infinity",               0 },
        { "Infinit",                0 },
        { "Infinityx",              0 },
        { "[true false]",           6 },
        // Strings
        { "\"abc",                  4 },
        { "\"",                     1 },
        { "\"\\",                   2 },
        { "\"a\\x\"",               2 },
        { "\"a\\'\"",               2 },
        { "\"a\\u12\"",             2 },
        { "\"a\\u12G4\"",           2 },
        { "\"a\\u",                 2 },
        { "\"\t\"",                 1 },
        { "\"a\nb\"",               2 },
        { "\"a\x01\"",              2 },
        { "'a'",                    0 },
        // Not JSON
        { "/* c */ 1",              0 },
        { "1 // c",                 2 },
        { "\xEF\xBB\xBF",           3 },
        { "\xEF\xBB\xBF\xEF\xBB\xBF" "1", 3 },
        { "\xEF\xBB" "1",           0 },
        { " \xEF\xBB\xBF" "1",      1 },
    };
    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        parsed_t p = parse(cases[i].text);
        EXPECT_TRUE_MSG(p.doc == NULL, cases[i].text);
        EXPECT_EQ_MSG(cases[i].offset, p.err.offset, cases[i].text);
        EXPECT_TRUE_MSG(p.err.message != NULL && p.err.message[0] != '\0', cases[i].text);
        release(&p);
    }

    // A NUL is not whitespace, and a NULL text is an empty one
    {
        const char text[] = "[1,\0 2]";
        parsed_t p = parse_n(text, sizeof(text) - 1);
        EXPECT_TRUE(p.doc == NULL);
        EXPECT_EQ(3, p.err.offset);
        release(&p);
    }
    {
        md_json_error_t err = {0};
        str_t null_text = { NULL, 0 };
        EXPECT_TRUE(md_json_parse(null_text, md_get_heap_allocator(), &err) == NULL);
        EXPECT_EQ(0, err.offset);
    }
}

UTEST(json, error_position) {
    const char* text =
        "{\n"
        "  \"a\": [1,\n"
        "    2,,\n"
        "  ]\n"
        "}\n";
    parsed_t p = parse(text);
    ASSERT_TRUE(p.doc == NULL);
    EXPECT_EQ(3, p.err.line);
    EXPECT_EQ(7, p.err.column);
    EXPECT_EQ((size_t)(strstr(text, ",,") - text) + 1, p.err.offset);
    release(&p);

    p = parse("[1,\r\n\"abc");
    ASSERT_TRUE(p.doc == NULL);
    EXPECT_EQ(2, p.err.line);
    EXPECT_EQ(5, p.err.column);
    release(&p);
}

// A view into a larger buffer is parsed as exactly that view
UTEST(json, reads_only_its_view) {
    static const struct {
        const char* buffer;
        size_t      len;
        bool        valid;
    } cases[] = {
        { "[1,2]3]",  5, true  },
        { "123",      2, true  },
        { "true",     3, false },
        { "\"ab\"",   3, false },
        { "1.5",      2, false },
        { "1e5",      2, false },
        { "-Infinity", 5, false },
        { "{\"a\":1}", 6, false },
        { "\"\\u0041\"", 6, false },
    };
    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        parsed_t p = parse_n(cases[i].buffer, cases[i].len);
        EXPECT_EQ_MSG(cases[i].valid, p.doc != NULL, cases[i].buffer);
        release(&p);
    }

    parsed_t p = parse_n("123", 2);
    int64_t v = 0;
    EXPECT_TRUE(md_json_i64(&v, md_json_root(p.doc)));
    EXPECT_EQ(12, v);
    release(&p);
}

// No recursion: the depth is limited by nothing but memory
UTEST(json, deep_nesting) {
    const size_t depth = 200000;
    char* text = malloc(depth * 2 + 1);
    ASSERT_TRUE(text);
    memset(text, '[', depth);
    memset(text + depth, ']', depth);
    text[2 * depth] = '\0';

    parsed_t p = parse(text);
    ASSERT_TRUE(p.doc != NULL);
    md_json_val_t v = md_json_root(p.doc);
    size_t levels = 0;
    while (md_json_count(v) == 1) {
        v = md_json_first(v);
        levels += 1;
    }
    EXPECT_EQ(depth - 1, levels);
    EXPECT_EQ(MD_JSON_TYPE_ARRAY, md_json_type(v));
    EXPECT_EQ(0, md_json_count(v));
    release(&p);

    text[2 * depth - 1] = '\0';
    p = parse(text);
    EXPECT_TRUE(p.doc == NULL);
    EXPECT_EQ(2 * depth - 1, p.err.offset);
    release(&p);
    free(text);
}

// The allocation is sized by counting separators, including any inside strings
UTEST(json, separators_inside_strings) {
    parsed_t p = parse("[\"a,b:c[d{e\", \",,,,\", {\"k:,[{\": [\"]}\"]}]");
    ASSERT_TRUE(p.doc != NULL);
    md_json_val_t root = md_json_root(p.doc);
    EXPECT_EQ(3, md_json_count(root));
    EXPECT_TRUE(md_json_string_eq(md_json_at(root, 0), str_of("a,b:c[d{e")));
    EXPECT_TRUE(md_json_string_eq(md_json_at(root, 1), str_of(",,,,")));
    EXPECT_TRUE(md_json_string_eq(md_json_at(get(md_json_at(root, 2), "k:,[{"), 0), str_of("]}")));
    release(&p);
}

// ---------------------------------------------------------------------------
// Navigation
// ---------------------------------------------------------------------------

static const char* water_text =
    "{\n"
    "  \"name\": \"water\",\n"
    "  \"atoms\": [\n"
    "    {\"element\": \"O\", \"xyz\": [0.0, 0.0, 0.1173]},\n"
    "    {\"element\": \"H\", \"xyz\": [0.0, 0.7572, -0.4692]},\n"
    "    {\"element\": \"H\", \"xyz\": [0.0, -0.7572, -0.4692]}\n"
    "  ],\n"
    "  \"empty_array\": [],\n"
    "  \"empty_object\": {},\n"
    "  \"nested\": {\"a\": {\"b\": {\"c\": [1, [2, [3]], {}]}}},\n"
    "  \"flag\": true,\n"
    "  \"nothing\": null,\n"
    "  \"dup\": 1,\n"
    "  \"dup\": [2]\n"
    "}\n";

UTEST(json, navigation) {
    parsed_t p = parse(water_text);
    ASSERT_TRUE(p.doc != NULL);
    const md_json_val_t root = md_json_root(p.doc);

    EXPECT_EQ(MD_JSON_TYPE_OBJECT, md_json_type(root));
    EXPECT_EQ(9, md_json_count(root));

    EXPECT_TRUE(md_json_string_eq(get(root, "name"), str_of("water")));

    const md_json_val_t atoms = get(root, "atoms");
    EXPECT_EQ(MD_JSON_TYPE_ARRAY, md_json_type(atoms));
    ASSERT_EQ(3, md_json_count(atoms));
    const char* elements[] = { "O", "H", "H" };
    const double y[] = { 0.0, 0.7572, -0.7572 };
    size_t i = 0;
    for (md_json_val_t a = md_json_first(atoms); md_json_valid(a); a = md_json_next(a), ++i) {
        ASSERT_LT(i, 3);
        EXPECT_EQ(md_json_at(atoms, i).idx, a.idx);
        EXPECT_TRUE(md_json_string_eq(get(a, "element"), str_of(elements[i])));
        double xyz[4] = { -1, -1, -1, -1 };
        EXPECT_EQ(3, md_json_extract_f64(xyz, 4, get(a, "xyz")));
        EXPECT_EQ(bits_of(y[i]), bits_of(xyz[1]));
        EXPECT_EQ(bits_of(-1.0), bits_of(xyz[3]));
        // An element of an array is no member, and has no key
        EXPECT_FALSE(md_json_valid(md_json_key(a)));
    }
    EXPECT_EQ(3, i);

    EXPECT_EQ(MD_JSON_TYPE_ARRAY,  md_json_type(get(root, "empty_array")));
    EXPECT_EQ(0, md_json_count(get(root, "empty_array")));
    EXPECT_FALSE(md_json_valid(md_json_first(get(root, "empty_array"))));
    EXPECT_EQ(MD_JSON_TYPE_OBJECT, md_json_type(get(root, "empty_object")));
    EXPECT_FALSE(md_json_valid(md_json_first(get(root, "empty_object"))));

    // Nested containers are stepped over whole
    const md_json_val_t c = get(get(get(get(root, "nested"), "a"), "b"), "c");
    ASSERT_EQ(3, md_json_count(c));
    int64_t v = 0;
    EXPECT_TRUE(md_json_i64(&v, md_json_at(c, 0)));
    EXPECT_EQ(1, v);
    EXPECT_TRUE(md_json_i64(&v, md_json_at(md_json_at(md_json_at(c, 1), 1), 0)));
    EXPECT_EQ(3, v);
    EXPECT_EQ(MD_JSON_TYPE_OBJECT, md_json_type(md_json_at(c, 2)));
    EXPECT_EQ(md_json_at(c, 2).idx, md_json_next(md_json_at(c, 1)).idx);
    EXPECT_FALSE(md_json_valid(md_json_next(md_json_at(c, 2))));
    EXPECT_FALSE(md_json_valid(md_json_at(c, 3)));

    bool flag = false;
    EXPECT_TRUE(md_json_bool(&flag, get(root, "flag")));
    EXPECT_TRUE(flag);
    EXPECT_EQ(MD_JSON_TYPE_NULL, md_json_type(get(root, "nothing")));
    EXPECT_FALSE(md_json_bool(&flag, get(root, "nothing")));

    // A repeated key: the last one is found, as in Python and JavaScript
    EXPECT_EQ(MD_JSON_TYPE_ARRAY, md_json_type(get(root, "dup")));

    // Members in document order, by walking and by index, each with its key
    const char* keys[] = { "name", "atoms", "empty_array", "empty_object", "nested", "flag", "nothing", "dup", "dup" };
    i = 0;
    for (md_json_val_t m = md_json_first(root); md_json_valid(m); m = md_json_next(m), ++i) {
        ASSERT_LT(i, ARRAY_SIZE(keys));
        EXPECT_EQ(md_json_at(root, i).idx, m.idx);
        const md_json_val_t k = md_json_key(m);
        EXPECT_EQ(MD_JSON_TYPE_STRING, md_json_type(k));
        EXPECT_TRUE(md_json_string_eq(k, str_of(keys[i])));
        // A key is not a member of its own
        EXPECT_FALSE(md_json_valid(md_json_next(k)));
        EXPECT_FALSE(md_json_valid(md_json_key(k)));
    }
    EXPECT_EQ(ARRAY_SIZE(keys), i);

    release(&p);
}

// NONE answers everything as an empty value would, so a chain of lookups needs one check
UTEST(json, none_propagates) {
    parsed_t p = parse(water_text);
    ASSERT_TRUE(p.doc != NULL);
    const md_json_val_t root = md_json_root(p.doc);

    const md_json_val_t nones[] = {
        { 0, 0 },
        md_json_root(NULL),
        get(root, "missing"),
        get(root, "Name"),
        get(get(root, "missing"), "deeper"),
        get(get(root, "name"), "x"),         // lookup in a string
        get(get(root, "atoms"), "0"),        // key lookup in an array
        md_json_at(get(root, "atoms"), 3),
        md_json_at(get(root, "atoms"), SIZE_MAX),
        md_json_at(get(root, "flag"), 0),
        md_json_first(get(root, "flag")),
        md_json_at(md_json_at(md_json_at(root, 100), 0), 0),
    };
    for (size_t i = 0; i < ARRAY_SIZE(nones); ++i) {
        const md_json_val_t v = nones[i];
        EXPECT_EQ(MD_JSON_TYPE_NONE, md_json_type(v));
        EXPECT_FALSE(md_json_valid(v));
        EXPECT_EQ(0, md_json_count(v));
        EXPECT_FALSE(md_json_valid(md_json_first(v)));
        EXPECT_FALSE(md_json_valid(md_json_next(v)));
        EXPECT_FALSE(md_json_valid(md_json_key(v)));
        EXPECT_FALSE(md_json_valid(md_json_get(v, str_of("a"))));
        EXPECT_FALSE(md_json_valid(md_json_at(v, 0)));

        double d = 42.0;
        int64_t n = 42;
        bool b = true;
        EXPECT_FALSE(md_json_f64(&d, v));
        EXPECT_FALSE(md_json_i64(&n, v));
        EXPECT_FALSE(md_json_bool(&b, v));
        EXPECT_EQ(bits_of(42.0), bits_of(d));
        EXPECT_EQ(42, n);
        EXPECT_TRUE(b);

        EXPECT_EQ(0, md_json_extract_f64(&d, 1, v));
        EXPECT_EQ(0, md_json_string_raw(v).len);
        EXPECT_FALSE(md_json_string_eq(v, str_of("")));
        char buf[4] = "xyz";
        EXPECT_EQ(0, md_json_string_copy(buf, sizeof(buf), v));
        EXPECT_EQ('\0', buf[0]);
        EXPECT_EQ(0, md_json_string(v, md_get_heap_allocator()).len);
    }

    release(&p);
}

// ---------------------------------------------------------------------------
// Values
// ---------------------------------------------------------------------------

UTEST(json, numbers) {
    static const char* texts[] = {
        "0", "-0", "1", "-1", "0.5", "-0.0e-0", "1e5", "1E5", "1e-5", "1e+5", "123456789012345678901234567890",
        "0.1", "-248.271233", "4.370019985270747", "1.396125e-30", "1e400", "-1e400", "1e-400", "2.2250738585072011e-308",
        "4.9406564584124654e-324", "1.7976931348623157e308", "NaN", "Infinity", "-Infinity",
    };
    for (size_t i = 0; i < ARRAY_SIZE(texts); ++i) {
        parsed_t p = parse(texts[i]);
        ASSERT_TRUE_MSG(p.doc != NULL, texts[i]);
        double expected = 0, got = 0;
        md_parse_f64(&expected, str_of(texts[i]));
        EXPECT_TRUE_MSG(md_json_f64(&got, md_json_root(p.doc)), texts[i]);
        EXPECT_EQ_MSG(bits_of(expected), bits_of(got), texts[i]);
        release(&p);
    }

    static const struct {
        const char* text;
        bool        ok;
        int64_t     value;
    } ints[] = {
        { "0",                      true,  0 },
        { "-0",                     true,  0 },
        { "42",                     true,  42 },
        { "-42",                    true,  -42 },
        { "9223372036854775807",    true,  INT64_MAX },
        { "-9223372036854775808",   true,  INT64_MIN },
        { "9223372036854775808",    false, 0 },
        { "-9223372036854775809",   false, 0 },
        { "99999999999999999999",   false, 0 },
        { "3.0",                    true,  3 },
        { "-3.0",                   true,  -3 },
        { "1e3",                    true,  1000 },
        { "1.5e1",                  true,  15 },
        { "-0.0",                   true,  0 },
        { "9007199254740991.0",     true,  9007199254740991 },
        { "9007199254740992.0",     false, 0 },
        { "9007199254740993.0",     false, 0 },
        { "2.5",                    false, 0 },
        { "1e-3",                   false, 0 },
        { "1e300",                  false, 0 },
        { "NaN",                    false, 0 },
        { "Infinity",               false, 0 },
        { "-Infinity",              false, 0 },
        { "\"1\"",                  false, 0 },
        { "true",                   false, 0 },
    };
    for (size_t i = 0; i < ARRAY_SIZE(ints); ++i) {
        parsed_t p = parse(ints[i].text);
        ASSERT_TRUE_MSG(p.doc != NULL, ints[i].text);
        int64_t v = 12345;
        EXPECT_EQ_MSG(ints[i].ok, md_json_i64(&v, md_json_root(p.doc)), ints[i].text);
        EXPECT_EQ_MSG(ints[i].ok ? ints[i].value : 12345, v, ints[i].text);
        release(&p);
    }

    // Booleans are only true and false
    parsed_t p = parse("[true, false, null, 0, 1, \"true\"]");
    ASSERT_TRUE(p.doc != NULL);
    bool b = false;
    EXPECT_TRUE(md_json_bool(&b, md_json_at(md_json_root(p.doc), 0)));
    EXPECT_TRUE(b);
    EXPECT_TRUE(md_json_bool(&b, md_json_at(md_json_root(p.doc), 1)));
    EXPECT_FALSE(b);
    for (size_t i = 2; i < 6; ++i) {
        b = true;
        EXPECT_FALSE(md_json_bool(&b, md_json_at(md_json_root(p.doc), i)));
        EXPECT_TRUE(b);
    }
    double d = 7.0;
    EXPECT_FALSE(md_json_f64(&d, md_json_at(md_json_root(p.doc), 0)));
    EXPECT_FALSE(md_json_f64(&d, md_json_at(md_json_root(p.doc), 5)));
    EXPECT_EQ(bits_of(7.0), bits_of(d));
    release(&p);
}

UTEST(json, extract) {
    parsed_t p = parse("{\"a\": [1, -2.5, 3e2, NaN, 5], \"b\": [1, \"x\", 3], \"c\": [], \"d\": {\"x\": 1}, \"e\": [[1], 2]}");
    ASSERT_TRUE(p.doc != NULL);
    const md_json_val_t root = md_json_root(p.doc);

    double d[8];
    float  f[8];
    EXPECT_EQ(5, md_json_extract_f64(d, 8, get(root, "a")));
    EXPECT_EQ(bits_of(-2.5), bits_of(d[1]));
    EXPECT_EQ(bits_of(300.0), bits_of(d[2]));
    EXPECT_TRUE((bits_of(d[3]) & 0x7FFFFFFFFFFFFFFFull) > 0x7FF0000000000000ull);
    EXPECT_EQ(5, md_json_extract_f32(f, 8, get(root, "a")));
    EXPECT_EQ(300.0f, f[2]);
    EXPECT_EQ(5.0f, f[4]);

    // Capped, and stopping at the first element that is no number
    EXPECT_EQ(2, md_json_extract_f64(d, 2, get(root, "a")));
    EXPECT_EQ(0, md_json_extract_f64(d, 0, get(root, "a")));
    EXPECT_EQ(0, md_json_extract_f64(NULL, 0, get(root, "a")));
    EXPECT_EQ(1, md_json_extract_f64(d, 8, get(root, "b")));
    EXPECT_EQ(1, md_json_extract_f32(f, 8, get(root, "b")));
    EXPECT_EQ(0, md_json_extract_f64(d, 8, get(root, "c")));
    EXPECT_EQ(0, md_json_extract_f64(d, 8, get(root, "d")));    // an object is no array
    EXPECT_EQ(0, md_json_extract_f64(d, 8, get(root, "e")));
    release(&p);
}

UTEST(json, strings) {
    parsed_t p = parse(
        "[\"plain\","
        " \"a\\\"b\\\\c\\/d\\b\\f\\n\\r\\te\","
        " \"\\u0041\\u00e9\\u00E9\\u20ac\\ud83d\\ude00\","
        " \"\\ud83d\","
        " \"\\ude00x\","
        " \"\\ud83d\\u0041\","
        " \"\\ud83d\\ud83d\\ude00\","
        " \"\\u0000\","
        " \"\xC3\xA9\xE2\x82\xAC\","
        " \"\"]");
    ASSERT_TRUE(p.doc != NULL);
    const md_json_val_t root = md_json_root(p.doc);
    ASSERT_EQ(10, md_json_count(root));

    static const struct {
        const char* decoded;
        size_t      len;
    } expected[] = {
        { "plain", 5 },
        { "a\"b\\c/d\b\f\n\r\te", 13 },
        { "A\xC3\xA9\xC3\xA9\xE2\x82\xAC\xF0\x9F\x98\x80", 12 },
        { "\xEF\xBF\xBD", 3 },                          // a high surrogate alone
        { "\xEF\xBF\xBDx", 4 },                         // a low surrogate alone
        { "\xEF\xBF\xBD" "A", 4 },                      // a high one followed by something else
        { "\xEF\xBF\xBD\xF0\x9F\x98\x80", 7 },          // two highs, then a low
        { "\0", 1 },
        { "\xC3\xA9\xE2\x82\xAC", 5 },                  // UTF-8 as it was written
        { "", 0 },
    };

    md_allocator_i* heap = md_get_heap_allocator();
    for (size_t i = 0; i < ARRAY_SIZE(expected); ++i) {
        const md_json_val_t v = md_json_at(root, i);
        const str_t want = { expected[i].decoded, expected[i].len };

        EXPECT_TRUE(md_json_string_eq(v, want));
        if (want.len) {
            const str_t shorter = { want.ptr, want.len - 1 };
            EXPECT_FALSE(md_json_string_eq(v, shorter));
        }
        const str_t longer = str_of("plain!");
        EXPECT_FALSE(md_json_string_eq(v, longer));

        char buf[32];
        memset(buf, 'x', sizeof(buf));
        const size_t n = md_json_string_copy(buf, sizeof(buf), v);
        EXPECT_EQ(want.len, n);
        EXPECT_TRUE(memcmp(buf, want.ptr, want.len) == 0);
        EXPECT_EQ('\0', buf[n]);

        const str_t s = md_json_string(v, heap);
        EXPECT_EQ(want.len, s.len);
        if (s.len) {
            EXPECT_TRUE(memcmp(s.ptr, want.ptr, want.len) == 0);
            EXPECT_EQ('\0', s.ptr[s.len]);
            str_free(s, heap);
        }
    }

    // Raw is the text between the quotes
    EXPECT_TRUE(str_eq(md_json_string_raw(md_json_at(root, 0)), str_of("plain")));
    EXPECT_TRUE(str_eq(md_json_string_raw(md_json_at(root, 3)), str_of("\\ud83d")));
    EXPECT_EQ(0, md_json_string_raw(md_json_at(root, 9)).len);
    release(&p);
}

// A copy that does not fit is cut between characters, both for UTF-8 written as is and decoded from escapes
UTEST(json, string_copy_cuts_between_characters) {
    const char* texts[] = {
        "\"a\xC3\xA9\xE2\x82\xAC\xF0\x9F\x98\x80\"",
        "\"a\\u00e9\\u20ac\\ud83d\\ude00\"",
    };
    // a (1) e-acute (2) euro (3) grinning face (4): what fits in cap - 1 bytes
    const size_t expected_len[] = { 0, 0, 1, 1, 3, 3, 3, 6, 6, 6, 6, 10, 10 };
    for (size_t t = 0; t < ARRAY_SIZE(texts); ++t) {
        parsed_t p = parse(texts[t]);
        ASSERT_TRUE(p.doc != NULL);
        for (size_t cap = 0; cap < ARRAY_SIZE(expected_len); ++cap) {
            char buf[16];
            memset(buf, 'x', sizeof(buf));
            const size_t n = md_json_string_copy(buf, cap, md_json_root(p.doc));
            EXPECT_EQ(expected_len[cap], n);
            if (cap) {
                EXPECT_EQ('\0', buf[n]);
            } else {
                EXPECT_EQ('x', buf[0]);
            }
        }
        release(&p);
    }
}

UTEST(json, keys_compare_decoded) {
    parsed_t p = parse("{\"a\\u0062c\": 1, \"\": 2, \"\\\"q\\\"\": 3, \"x\\/y\": 4}");
    ASSERT_TRUE(p.doc != NULL);
    const md_json_val_t root = md_json_root(p.doc);
    int64_t v = 0;
    EXPECT_TRUE(md_json_i64(&v, get(root, "abc")));
    EXPECT_EQ(1, v);
    EXPECT_FALSE(md_json_valid(get(root, "a\\u0062c")));
    EXPECT_FALSE(md_json_valid(get(root, "ab")));
    EXPECT_FALSE(md_json_valid(get(root, "abcd")));
    EXPECT_TRUE(md_json_i64(&v, get(root, "")));
    EXPECT_EQ(2, v);
    EXPECT_TRUE(md_json_i64(&v, get(root, "\"q\"")));
    EXPECT_EQ(3, v);
    EXPECT_TRUE(md_json_i64(&v, get(root, "x/y")));
    EXPECT_EQ(4, v);
    release(&p);
}

// ---------------------------------------------------------------------------
// Random documents: generated together with what a walk over them must find
// ---------------------------------------------------------------------------

typedef struct text_buf_t {
    char*  ptr;
    size_t len;
    size_t cap;
} text_buf_t;

static void buf_push(text_buf_t* b, const char* s, size_t n) {
    if (b->len + n + 1 > b->cap) {
        b->cap = (b->len + n + 1) * 2;
        b->ptr = realloc(b->ptr, b->cap);
    }
    memcpy(b->ptr + b->len, s, n);
    b->len += n;
    b->ptr[b->len] = '\0';
}

static void buf_pushc(text_buf_t* b, const char* s) {
    buf_push(b, s, strlen(s));
}

// Pre-order: a value is an entry, a member's key is the entry in front of its value
typedef struct expect_t {
    md_json_type_t type;
    bool     is_key;
    bool     bool_value;
    size_t   count;
    size_t   num_off, num_len;      // the number's text, in the generated document
    size_t   str_off, str_len;      // the decoded string, in the strings buffer
} expect_t;

typedef struct gen_t {
    uint64_t    rng;
    text_buf_t  text;
    text_buf_t  strings;
    expect_t*   exp;
    size_t      num_exp, cap_exp;
} gen_t;

static size_t gen_expect(gen_t* g, expect_t e) {
    if (g->num_exp == g->cap_exp) {
        g->cap_exp = g->cap_exp ? g->cap_exp * 2 : 256;
        g->exp = realloc(g->exp, g->cap_exp * sizeof(expect_t));
    }
    g->exp[g->num_exp] = e;
    return g->num_exp++;
}

static void gen_ws(gen_t* g) {
    static const char* ws[] = { "", "", "", " ", "\n", "  ", "\t", "\r\n    " };
    buf_pushc(&g->text, ws[rng_range(&g->rng, ARRAY_SIZE(ws))]);
}

static size_t utf8_of(char out[4], uint32_t cp) {
    if (cp < 0x80)    { out[0] = (char)cp; return 1; }
    if (cp < 0x800)   { out[0] = (char)(0xC0 | (cp >> 6)); out[1] = (char)(0x80 | (cp & 0x3F)); return 2; }
    if (cp < 0x10000) { out[0] = (char)(0xE0 | (cp >> 12)); out[1] = (char)(0x80 | ((cp >> 6) & 0x3F)); out[2] = (char)(0x80 | (cp & 0x3F)); return 3; }
    out[0] = (char)(0xF0 | (cp >> 18)); out[1] = (char)(0x80 | ((cp >> 12) & 0x3F)); out[2] = (char)(0x80 | ((cp >> 6) & 0x3F)); out[3] = (char)(0x80 | (cp & 0x3F));
    return 4;
}

static void gen_string(gen_t* g, bool is_key) {
    expect_t e = { MD_JSON_TYPE_STRING, is_key };
    e.str_off = g->strings.len;
    buf_pushc(&g->text, "\"");
    const uint32_t len = rng_range(&g->rng, 4) == 0 ? rng_range(&g->rng, 40) : rng_range(&g->rng, 6);
    for (uint32_t i = 0; i < len; ++i) {
        char enc[16], dec[4];
        size_t enc_len, dec_len;
        uint32_t cp;
        switch (rng_range(&g->rng, 8)) {
        case 0:  cp = rng_range(&g->rng, 0x20); break;                            // control: escaped
        case 1:  cp = "\"\\/"[rng_range(&g->rng, 3)]; break;
        case 2:  cp = 0x80 + rng_range(&g->rng, 0x780); break;                    // two bytes
        case 3:  cp = 0x800 + rng_range(&g->rng, 0xD000 - 0x800); break;         // three, below the surrogates
        case 4:  cp = 0x10000 + rng_range(&g->rng, 0x100000); break;              // four
        default: cp = 0x20 + rng_range(&g->rng, 0x5F); break;                     // printable ASCII
        }
        dec_len = utf8_of(dec, cp);
        const bool escape = cp < 0x20 || cp == '"' || cp == '\\' || rng_range(&g->rng, 3) == 0;
        if (!escape) {
            memcpy(enc, dec, dec_len);
            enc_len = dec_len;
        } else if (cp == '"' || cp == '\\' || cp == '/' || cp == '\b' || cp == '\f' || cp == '\n' || cp == '\r' || cp == '\t') {
            const char* short_esc = cp == '"' ? "\\\"" : cp == '\\' ? "\\\\" : cp == '/' ? "\\/" : cp == '\b' ? "\\b" :
                                    cp == '\f' ? "\\f" : cp == '\n' ? "\\n" : cp == '\r' ? "\\r" : "\\t";
            enc_len = (size_t)snprintf(enc, sizeof(enc), "%s", short_esc);
        } else if (cp < 0x10000) {
            enc_len = (size_t)snprintf(enc, sizeof(enc), rng_range(&g->rng, 2) ? "\\u%04x" : "\\u%04X", cp);
        } else {
            const uint32_t v = cp - 0x10000;
            enc_len = (size_t)snprintf(enc, sizeof(enc), "\\u%04x\\u%04X", 0xD800 + (v >> 10), 0xDC00 + (v & 0x3FF));
        }
        buf_push(&g->text, enc, enc_len);
        buf_push(&g->strings, dec, dec_len);
    }
    buf_pushc(&g->text, "\"");
    e.str_len = g->strings.len - e.str_off;
    gen_expect(g, e);
}

static void gen_number(gen_t* g) {
    expect_t e = { MD_JSON_TYPE_NUMBER };
    e.num_off = g->text.len;
    char num[64];
    int n = 0;
    switch (rng_range(&g->rng, 6)) {
    case 0:
        n = snprintf(num, sizeof(num), "%s", (const char*[]){ "NaN", "Infinity", "-Infinity", "0", "-0", "1e400" }[rng_range(&g->rng, 6)]);
        break;
    default: {
        if (rng_range(&g->rng, 2)) num[n++] = '-';
        const uint32_t int_digits = 1 + rng_range(&g->rng, 12);
        num[n++] = (char)('1' + rng_range(&g->rng, 9));
        for (uint32_t i = 1; i < int_digits; ++i) num[n++] = (char)('0' + rng_range(&g->rng, 10));
        if (rng_range(&g->rng, 2)) {
            num[n++] = '.';
            const uint32_t frac_digits = 1 + rng_range(&g->rng, 17);
            for (uint32_t i = 0; i < frac_digits; ++i) num[n++] = (char)('0' + rng_range(&g->rng, 10));
        }
        if (rng_range(&g->rng, 3) == 0) {
            num[n++] = "eE"[rng_range(&g->rng, 2)];
            const uint32_t sign = rng_range(&g->rng, 3);
            if (sign) num[n++] = sign == 1 ? '+' : '-';
            n += snprintf(num + n, sizeof(num) - n, "%u", rng_range(&g->rng, 330));
        }
        break;
    }
    }
    buf_push(&g->text, num, (size_t)n);
    e.num_len = (size_t)n;
    gen_expect(g, e);
}

static void gen_value(gen_t* g, int depth) {
    const uint32_t kind = depth > 6 ? 2 + rng_range(&g->rng, 4) : rng_range(&g->rng, 8);
    switch (kind) {
    case 0:
    case 1: {
        const bool object = kind == 1;
        const size_t self = gen_expect(g, (expect_t){ object ? MD_JSON_TYPE_OBJECT : MD_JSON_TYPE_ARRAY });
        const uint32_t count = rng_range(&g->rng, 4) == 0 ? 0 : 1 + rng_range(&g->rng, 6);
        buf_pushc(&g->text, object ? "{" : "[");
        for (uint32_t i = 0; i < count; ++i) {
            if (i) buf_pushc(&g->text, ",");
            gen_ws(g);
            if (object) {
                gen_string(g, true);
                gen_ws(g);
                buf_pushc(&g->text, ":");
                gen_ws(g);
            }
            gen_value(g, depth + 1);
            gen_ws(g);
        }
        buf_pushc(&g->text, object ? "}" : "]");
        g->exp[self].count = count;
        break;
    }
    case 2: gen_string(g, false); break;
    case 3:
    case 4: gen_number(g); break;
    default: {
        const uint32_t lit = rng_range(&g->rng, 3);
        buf_pushc(&g->text, lit == 0 ? "true" : lit == 1 ? "false" : "null");
        expect_t e = { lit == 2 ? MD_JSON_TYPE_NULL : MD_JSON_TYPE_BOOL };
        e.bool_value = lit == 0;
        gen_expect(g, e);
        break;
    }
    }
}

typedef struct check_t {
    const gen_t* g;
    size_t       next;      // the next expectation
    size_t       failures;
} check_t;

static bool check_string(const gen_t* g, const expect_t* e, md_json_val_t v) {
    const str_t want = { g->strings.ptr + e->str_off, e->str_len };
    if (!md_json_string_eq(v, want)) return false;
    char buf[512];
    const size_t n = md_json_string_copy(buf, sizeof(buf), v);
    if (n != want.len || (n && memcmp(buf, want.ptr, n) != 0)) return false;
    const str_t s = md_json_string(v, md_get_heap_allocator());
    const bool ok = s.len == want.len && (s.len == 0 || memcmp(s.ptr, want.ptr, s.len) == 0);
    if (s.len) str_free(s, md_get_heap_allocator());
    return ok;
}

static void check_value(check_t* c, md_json_val_t v) {
    const gen_t* g = c->g;
    if (c->next >= g->num_exp) { c->failures++; return; }
    const expect_t* e = &g->exp[c->next++];
    if (md_json_type(v) != e->type) { c->failures++; return; }

    switch (e->type) {
    case MD_JSON_TYPE_STRING:
        if (!check_string(g, e, v)) c->failures++;
        break;
    case MD_JSON_TYPE_NUMBER: {
        const str_t text = { g->text.ptr + e->num_off, e->num_len };
        double want = 0, got = 0;
        md_parse_f64(&want, text);
        if (!md_json_f64(&got, v) || bits_of(want) != bits_of(got)) c->failures++;
        if (md_json_string_raw(v).len != 0) c->failures++;
        break;
    }
    case MD_JSON_TYPE_BOOL: {
        bool b = !e->bool_value;
        if (!md_json_bool(&b, v) || b != e->bool_value) c->failures++;
        break;
    }
    case MD_JSON_TYPE_ARRAY:
    case MD_JSON_TYPE_OBJECT: {
        const bool object = e->type == MD_JSON_TYPE_OBJECT;
        if (md_json_count(v) != e->count) { c->failures++; return; }
        size_t i = 0;
        for (md_json_val_t m = md_json_first(v); md_json_valid(m); m = md_json_next(m), ++i) {
            if (i >= e->count) { c->failures++; return; }
            if (md_json_at(v, i).idx != m.idx) c->failures++;
            if (object) {
                const expect_t* k = &g->exp[c->next++];
                if (!k->is_key || !check_string(g, k, md_json_key(m))) c->failures++;
                // A lookup by this key finds this member, unless a later one repeats it
                const str_t key = { g->strings.ptr + k->str_off, k->str_len };
                const md_json_val_t found = md_json_get(v, key);
                bool repeated = false;
                for (md_json_val_t later = md_json_next(m); md_json_valid(later); later = md_json_next(later)) {
                    if (md_json_string_eq(md_json_key(later), key)) repeated = true;
                }
                if (!repeated && found.idx != m.idx) c->failures++;
            } else if (md_json_valid(md_json_key(m))) {
                c->failures++;
            }
            check_value(c, m);
        }
        if (i != e->count) c->failures++;
        break;
    }
    default:
        break;
    }
}

UTEST(json, random_documents) {
    gen_t g = {0};
    g.rng = 0x9E3779B97F4A7C15ull;
    size_t total_bytes = 0;
    for (int doc = 0; doc < 2000; ++doc) {
        g.text.len = 0;
        g.strings.len = 0;
        g.num_exp = 0;
        buf_pushc(&g.text, "");
        buf_pushc(&g.strings, "");
        gen_ws(&g);
        gen_value(&g, 0);
        gen_ws(&g);
        total_bytes += g.text.len;

        parsed_t p = parse_n(g.text.ptr, g.text.len);
        if (!p.doc) {
            printf("  failed to parse: %s at %zu in\n%s\n", p.err.message, p.err.offset, g.text.ptr);
            ASSERT_TRUE(false);
        }
        check_t c = { &g, 0, 0 };
        check_value(&c, md_json_root(p.doc));
        if (c.failures || c.next != g.num_exp) {
            printf("  %zu mismatches in\n%s\n", c.failures, g.text.ptr);
        }
        ASSERT_EQ(0, c.failures);
        ASSERT_EQ(g.num_exp, c.next);
        release(&p);
    }
    EXPECT_GT(total_bytes, 100000);
    free(g.text.ptr);
    free(g.strings.ptr);
    free(g.exp);
}

// Every way of breaking a document is either an error or a document that can be walked in full
static void walk_all(md_json_val_t v, size_t* visited) {
    *visited += 1;
    double d;
    int64_t i;
    bool b;
    char buf[8];
    md_json_f64(&d, v);
    md_json_i64(&i, v);
    md_json_bool(&b, v);
    md_json_string_copy(buf, sizeof(buf), v);
    md_json_string_eq(v, str_of("abc"));
    md_json_get(v, str_of("a"));
    md_json_at(v, 2);
    md_json_extract_f64(&d, 1, v);
    md_json_key(v);
    for (md_json_val_t m = md_json_first(v); md_json_valid(m); m = md_json_next(m)) {
        md_json_string_eq(md_json_key(m), str_of("k"));
        walk_all(m, visited);
    }
}

UTEST(json, mutations) {
    gen_t g = {0};
    g.rng = 0xD1B54A32D192ED03ull;
    char* mutated = NULL;
    size_t accepted = 0, rejected = 0;
    for (int doc = 0; doc < 300; ++doc) {
        g.text.len = 0;
        g.strings.len = 0;
        g.num_exp = 0;
        buf_pushc(&g.text, "");
        buf_pushc(&g.strings, "");
        gen_value(&g, 0);
        mutated = realloc(mutated, g.text.len + 1);

        for (int m = 0; m < 40; ++m) {
            memcpy(mutated, g.text.ptr, g.text.len);
            size_t len = g.text.len;
            const uint32_t edits = 1 + rng_range(&g.rng, 3);
            for (uint32_t k = 0; k < edits && len; ++k) {
                const size_t at = rng_range(&g.rng, (uint32_t)len);
                switch (rng_range(&g.rng, 4)) {
                case 0: mutated[at] = "[]{}\",:\\0-eE.tfnu \x01\xC3"[rng_range(&g.rng, 21)]; break;
                case 1: mutated[at] = (char)rng_range(&g.rng, 256); break;
                case 2: len = at; break;                                                            // truncate
                case 3: memmove(mutated + at, mutated + at + 1, len - at - 1); len -= 1; break;     // delete
                }
            }
            parsed_t p = parse_n(mutated, len);
            if (p.doc) {
                size_t visited = 0;
                walk_all(md_json_root(p.doc), &visited);
                EXPECT_GT(visited, 0);
                accepted += 1;
            } else {
                EXPECT_LE(p.err.offset, len);
                EXPECT_TRUE(p.err.message != NULL);
                rejected += 1;
            }
            release(&p);
        }
    }
    EXPECT_GT(accepted, 0);
    EXPECT_GT(rejected, 0);
    free(mutated);
    free(g.text.ptr);
    free(g.strings.ptr);
    free(g.exp);
}

// ---------------------------------------------------------------------------
// Memory
// ---------------------------------------------------------------------------

typedef struct counting_t {
    md_allocator_i  iface;
    md_allocator_i* backing;
    int64_t         outstanding;    // bytes
    int64_t         allocations;
    bool            fail;
} counting_t;

static void* counting_realloc(struct md_allocator_o* inst, void* ptr, size_t old_size, size_t new_size, const char* file, size_t line) {
    counting_t* c = (counting_t*)inst;
    if (new_size && c->fail) return NULL;
    if (!ptr && new_size) c->allocations += 1;
    c->outstanding += (int64_t)new_size - (int64_t)(ptr ? old_size : 0);
    return c->backing->realloc(c->backing->inst, ptr, old_size, new_size, file, line);
}

static md_allocator_i* counting_init(counting_t* c) {
    memset(c, 0, sizeof(*c));
    c->backing = md_get_heap_allocator();
    c->iface.inst = (struct md_allocator_o*)c;
    c->iface.realloc = counting_realloc;
    return &c->iface;
}

UTEST(json, one_allocation) {
    counting_t c;
    md_allocator_i* alloc = counting_init(&c);

    md_json_t* doc = md_json_parse(str_of(water_text), alloc, NULL);
    ASSERT_TRUE(doc != NULL);
    EXPECT_EQ(1, c.allocations);
    // Bounded by 16 bytes per separator and a small constant
    size_t separators = 0;
    for (const char* s = water_text; *s; ++s) separators += (*s == ',' || *s == ':' || *s == '[' || *s == '{');
    EXPECT_LE(c.outstanding, (int64_t)(16 * separators + 128));

    // Reading allocates nothing
    md_json_string_copy((char[8]){0}, 8, get(md_json_root(doc), "name"));
    md_json_at(get(md_json_root(doc), "atoms"), 2);
    EXPECT_EQ(1, c.allocations);

    md_json_free(doc, alloc);
    EXPECT_EQ(0, c.outstanding);

    // A failed parse gives back what it took
    md_json_error_t err = {0};
    EXPECT_TRUE(md_json_parse(str_of("{\"a\": [1, 2,]}"), alloc, &err) == NULL);
    EXPECT_EQ(0, c.outstanding);

    // and an allocation that fails is an error, not a crash
    c.fail = true;
    EXPECT_TRUE(md_json_parse(str_of("[1, 2]"), alloc, &err) == NULL);
    EXPECT_TRUE(strstr(err.message, "memory") != NULL);
    EXPECT_EQ(0, c.outstanding);
}

UTEST(json, temp_arena) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp);
    md_json_t* doc = md_json_parse(str_of(water_text), arena, NULL);
    ASSERT_TRUE(doc != NULL);
    const str_t name = md_json_string(get(md_json_root(doc), "name"), arena);
    EXPECT_TRUE(str_eq(name, str_of("water")));
    md_temp_end(temp);
}

// ---------------------------------------------------------------------------
// Throughput, on a document shaped as the polarizable embedding VeloxChem writes
// ---------------------------------------------------------------------------

UTEST(json, perf) {
    uint64_t rng = 12345;
    text_buf_t b = {0};
    buf_pushc(&b, "{\n    \"classical_subsystems\": [\n        {\n            \"classical_fragments\": [\n");
    const int num_frag = 5000;
    for (int f = 0; f < num_frag; ++f) {
        char tmp[2048];
        snprintf(tmp, sizeof(tmp), "%s                {\n                    \"index\": %d,\n                    \"name\": \"HOH_pe\",\n                    \"atoms\": [\n", f ? ",\n" : "", f + 1);
        buf_pushc(&b, tmp);
        for (int a = 0; a < 3; ++a) {
            const double x = (double)(rng_next(&rng) >> 11) * (1.0 / 9007199254740992.0) * 40.0 - 20.0;
            const double y = (double)(rng_next(&rng) >> 11) * (1.0 / 9007199254740992.0) * 40.0 - 20.0;
            const double z = (double)(rng_next(&rng) >> 11) * (1.0 / 9007199254740992.0) * 40.0 - 20.0;
            snprintf(tmp, sizeof(tmp),
                "%s                        {\n"
                "                            \"index\": %d,\n"
                "                            \"element\": \"%s\",\n"
                "                            \"coordinate\": [\n                                %.17g,\n                                %.17g,\n                                %.17g\n                            ],\n"
                "                            \"multipoles\": {\"elements\": [%s]},\n"
                "                            \"polarizabilities\": {\"elements\": [0.0, 0.0, 0.0, 0.0, 5.73935, 0.0, 0.0, 5.73935, 0.0, 5.73935], \"order\": [1, 1]}\n"
                "                        }",
                a ? ",\n" : "", 3 * f + a + 1, a ? "H" : "O", x, y, z, a ? "0.33722" : "-0.67444");
            buf_pushc(&b, tmp);
        }
        buf_pushc(&b, "\n                    ]\n                }");
    }
    buf_pushc(&b, "\n            ]\n        }\n    ]\n}\n");

    md_tick_t t0 = md_tick_now();
    md_json_t* doc = md_json_parse((str_t){ b.ptr, b.len }, md_get_heap_allocator(), NULL);
    md_tick_t t1 = md_tick_now();
    ASSERT_TRUE(doc != NULL);

    double sum = 0;
    size_t num_atoms = 0;
    const md_json_val_t frags = get(md_json_at(get(md_json_root(doc), "classical_subsystems"), 0), "classical_fragments");
    for (md_json_val_t f = md_json_first(frags); md_json_valid(f); f = md_json_next(f)) {
        for (md_json_val_t a = md_json_first(get(f, "atoms")); md_json_valid(a); a = md_json_next(a)) {
            double xyz[3], q = 0, pol[10];
            md_json_extract_f64(xyz, 3, get(a, "coordinate"));
            md_json_f64(&q, md_json_at(get(get(a, "multipoles"), "elements"), 0));
            md_json_extract_f64(pol, 10, get(get(a, "polarizabilities"), "elements"));
            sum += xyz[0] + xyz[1] + xyz[2] + q + pol[4];
            num_atoms += 1;
        }
    }
    md_tick_t t2 = md_tick_now();
    EXPECT_EQ(3 * num_frag, num_atoms);

    const double mb = b.len / 1e6;
    const double parse_ms = md_tick_to_milliseconds(t1 - t0);
    const double walk_ms  = md_tick_to_milliseconds(t2 - t1);
    printf("  %.1f MB, %u values: parse %.2f ms (%.0f MB/s), walk and read numbers %.2f ms (sum %g)\n",
        mb, (unsigned)(num_atoms * 21), parse_ms, mb / (parse_ms * 1e-3), walk_ms, sum);

    md_json_free(doc, md_get_heap_allocator());
    free(b.ptr);
}

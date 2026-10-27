// The built-in fuzzer of the VIAL front end: every unittest run mutates the reference's examples, the golden inputs
// and the prelude, and checks every result against vial_fuzz_check.h. The mutations and their order come from a
// fixed seed, so a failure repeats; it prints the input that broke, as a C string.
//
//   VIAL_FUZZ_ITERATIONS=n    how many mutated inputs to check (default 3000)
//   VIAL_FUZZ_SEED=n          another sequence of mutations
//
// For long runs, use the libFuzzer target fuzz_vial.c, which runs the same check.

#include "utest.h"

#include <core/md_arena_allocator.h>
#include <core/md_common.h>
#include <core/md_str.h>

#include <stdio.h>
#include <stdlib.h>

#include "vial_spec.h"
#include "vial_fuzz_check.h"

#define VIAL_FUZZ_REFERENCE MD_UNITTEST_SOURCE_DIR "/../docs/vial_reference.md"
#define VIAL_FUZZ_PRELUDE   MD_UNITTEST_SOURCE_DIR "/../src/vial/prelude.vial"
#define VIAL_FUZZ_DICT      MD_UNITTEST_SOURCE_DIR "/vial_fuzz.dict"
#define VIAL_FUZZ_GOLDEN    MD_UNITTEST_SOURCE_DIR "/vial_golden/"

typedef struct fuzz_rng_t { uint64_t s; } fuzz_rng_t;

static uint64_t fuzz_next(fuzz_rng_t* r) {
    // xorshift64*
    r->s ^= r->s >> 12; r->s ^= r->s << 25; r->s ^= r->s >> 27;
    return r->s * 2685821657736338717ull;
}
static size_t fuzz_below(fuzz_rng_t* r, size_t n) { return n ? (size_t)(fuzz_next(r) % n) : 0; }

static uint64_t fuzz_env(const char* name, uint64_t def) {
    const char* v = getenv(name);
    return (v && *v) ? strtoull(v, NULL, 10) : def;
}

// The quoted entries of a libFuzzer dictionary, with \\, \" and \xNN unescaped
static str_t* fuzz_read_dict(str_t text, md_allocator_i* alloc) {
    str_t* words = NULL;
    str_t line;
    while (str_extract_line(&line, &text)) {
        line = str_trim(line);
        if (line.len < 2 || line.ptr[0] != '"' || line.ptr[line.len - 1] != '"') continue;
        char* w = md_alloc(alloc, line.len);
        size_t n = 0;
        for (size_t i = 1; i + 1 < line.len; ++i) {
            char c = line.ptr[i];
            if (c == '\\' && i + 2 < line.len) {
                const char e = line.ptr[++i];
                if (e == 'x' && i + 2 < line.len) {
                    char hex[3] = {line.ptr[i + 1], line.ptr[i + 2], 0};
                    c = (char)strtol(hex, NULL, 16);
                    i += 2;
                } else {
                    c = e;
                }
            }
            w[n++] = c;
        }
        md_array_push(words, ((str_t){w, n}), alloc);
    }
    return words;
}

typedef struct fuzz_buf_t { char* ptr; size_t len; } fuzz_buf_t;

static void fuzz_insert(fuzz_buf_t* b, size_t at, const char* src, size_t n, md_allocator_i* alloc) {
    char* p = md_alloc(alloc, b->len + n);
    memcpy(p, b->ptr, at);
    memcpy(p + at, src, n);
    memcpy(p + at + n, b->ptr + at, b->len - at);
    b->ptr = p;
    b->len += n;
}

static void fuzz_mutate(fuzz_buf_t* b, const str_t* seeds, const str_t* dict, fuzz_rng_t* r, md_allocator_i* alloc) {
    const size_t num_ops = 1 + fuzz_below(r, 6);
    for (size_t k = 0; k < num_ops; ++k) {
        const size_t at = fuzz_below(r, b->len + 1);
        switch (fuzz_below(r, 8)) {
        case 0:   // a word of the dictionary
        case 1: {
            const str_t w = dict[fuzz_below(r, md_array_size(dict))];
            fuzz_insert(b, at, w.ptr, w.len, alloc);
            break;
        }
        case 2:   // a random byte
            if (at < b->len) b->ptr[at] = (char)fuzz_next(r);
            break;
        case 3: { // a byte flipped
            if (at < b->len) b->ptr[at] ^= (char)(1u << fuzz_below(r, 8));
            break;
        }
        case 4: { // a range deleted
            const size_t n = fuzz_below(r, 16) + 1;
            const size_t end = at + n < b->len ? at + n : b->len;
            memmove(b->ptr + at, b->ptr + end, b->len - end);
            b->len -= end - at;
            break;
        }
        case 5: { // a range duplicated
            if (at >= b->len) break;
            const size_t n = fuzz_below(r, b->len - at < 64 ? b->len - at : 64) + 1;
            char tmp[64];
            memcpy(tmp, b->ptr + at, n);
            fuzz_insert(b, at, tmp, n, alloc);
            break;
        }
        case 6: { // a piece of another seed
            const str_t s = seeds[fuzz_below(r, md_array_size(seeds))];
            if (s.len == 0) break;
            const size_t beg = fuzz_below(r, s.len);
            const size_t n = fuzz_below(r, s.len - beg < 128 ? s.len - beg : 128) + 1;
            fuzz_insert(b, at, s.ptr + beg, n, alloc);
            break;
        }
        case 7: { // a line indented further, by spaces or a tab
            size_t ls = at;
            while (ls > 0 && b->ptr[ls - 1] != '\n') ls -= 1;
            if (fuzz_below(r, 4) == 0) fuzz_insert(b, ls, "\t", 1, alloc);
            else fuzz_insert(b, ls, "        ", 1 + fuzz_below(r, 8), alloc);
            break;
        }
        }
    }
}

static void fuzz_print_input(const char* ptr, size_t len) {
    printf("    \"");
    for (size_t i = 0; i < len; ++i) {
        const unsigned char c = (unsigned char)ptr[i];
        if (c == '\n') printf("\\n");
        else if (c == '"' || c == '\\') printf("\\%c", c);
        else if (c >= 32 && c < 127) printf("%c", c);
        else printf("\\x%02x\"\"", c);
    }
    printf("\"\n");
}

// Checks one input from an exact-size heap copy, so a read past its end is a sanitizer error, not a NUL
static const char* fuzz_check_copy(const char* ptr, size_t len, md_allocator_i* scratch) {
    char* copy = malloc(len ? len : 1);
    if (len) memcpy(copy, ptr, len);
    const char* err = vial_fuzz_check((str_t){copy, len}, scratch);
    free(copy);
    return err;
}

struct vial_fuzz {
    md_allocator_i* alloc;
    str_t* seeds;
    str_t* dict;
};

UTEST_F_SETUP(vial_fuzz) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    utest_fixture->alloc = alloc;
    utest_fixture->seeds = NULL;

    vial_spec_t spec = {0};
    str_t ref = load_textfile(STR_LIT(VIAL_FUZZ_REFERENCE), alloc);
    ASSERT_FALSE(str_empty(ref));
    vial_spec_read(&spec, ref, alloc);
    for (size_t i = 0; i < md_array_size(spec.blocks); ++i) md_array_push(utest_fixture->seeds, spec.blocks[i].src, alloc);

    const char* files[] = {VIAL_FUZZ_PRELUDE, VIAL_FUZZ_GOLDEN "expressions.vial", VIAL_FUZZ_GOLDEN "statements.vial",
                           VIAL_FUZZ_GOLDEN "refused.vial", VIAL_FUZZ_GOLDEN "errors.vial"};
    for (size_t i = 0; i < ARRAY_SIZE(files); ++i) {
        str_t s = load_textfile(str_from_cstr(files[i]), alloc);
        ASSERT_FALSE(str_empty(s));
        md_array_push(utest_fixture->seeds, s, alloc);
    }
    str_t dict = load_textfile(STR_LIT(VIAL_FUZZ_DICT), alloc);
    utest_fixture->dict = fuzz_read_dict(dict, alloc);
    ASSERT_GT(md_array_size(utest_fixture->dict), 50);
}

UTEST_F_TEARDOWN(vial_fuzz) {
    md_arena_allocator_destroy(utest_fixture->alloc);
}

// Every seed as it is, then every prefix of the golden error input: text cut off anywhere
UTEST_F(vial_fuzz, seeds_and_prefixes) {
    md_allocator_i* scratch = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    const str_t* seeds = utest_fixture->seeds;
    for (size_t i = 0; i < md_array_size(seeds); ++i) {
        const char* err = fuzz_check_copy(seeds[i].ptr, seeds[i].len, scratch);
        if (err) { printf("seed %zu: %s\n", i, err); fuzz_print_input(seeds[i].ptr, seeds[i].len); }
        EXPECT_TRUE(err == NULL);
        md_arena_allocator_reset(scratch);
    }
    for (size_t i = 0; i < md_array_size(seeds); ++i) {
        for (size_t n = 0; n <= seeds[i].len; n += 1 + n / 64) {
            const char* err = fuzz_check_copy(seeds[i].ptr, n, scratch);
            if (err) { printf("seed %zu cut at %zu: %s\n", i, n, err); fuzz_print_input(seeds[i].ptr, n); ASSERT_TRUE(err == NULL); }
            md_arena_allocator_reset(scratch);
        }
    }
    md_arena_allocator_destroy(scratch);
}

// Text no one writes by hand: deep nesting, long chains, many errors. Each must be reported (E0104 where it nests
// too deep) and stay within the stack and time a script gets.
UTEST_F(vial_fuzz, extremes) {
    md_allocator_i* scratch = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    typedef struct { const char* name; const char* head; const char* unit; const char* tail; int count; bool deep; char close; } extreme_t;
    const extreme_t cases[] = {
        {"parens",          "x = ",  "(",              "",        5000, true},
        {"parens closed",   "x = ",  "(",              "1",       5000, true, ')'},
        {"lists",           "x = ",  "[",              "",        5000, true},
        {"list comps",      "x = ",  "[a for a in ",   "",        2000, true},
        {"unary minus",     "x = ",  "-",              "1",       5000, true},
        {"not",             "x = ",  "not ",           "1",       5000, true},
        {"power",           "x = ",  "a**",            "1",       5000, true},
        {"conditional",     "x = ",  "a if b else ",   "c",       5000, true},
        {"lambda",          "x = ",  "lambda: ",       "1",       5000, true},
        {"starred",         "f(",    "*",              "a)",      5000, true},
        {"one-line ifs",    "",      "if a: ",         "pass",    5000, true},
        {"sum",             "x = a", " + a",           "",        5000, true},
        {"attributes",      "x = a", ".b",             "",        5000, true},
        {"calls",           "x = f", "()",             "",        5000, true},
        {"subscripts",      "x = a", "[1]",            "",        5000, true},
        {"comparisons",     "x = a", " < a",           "",        5000, false},
        {"tuple",           "x = ",  "1, ",            "",        5000, false},
        {"statements",      "",      "x = 1;",         "",        5000, false},
        {"error lines",     "",      "=\n",            "",        5000, false},
        {"lexer errors",    "",      "$ = (\n",        "",        5000, false},
        {"blank lines",     "",      "\n",             "",        5000, false},
        {"v1 calls",        "x = ",  "residue(1:4)",   "",        2000, false},
    };
    for (size_t c = 0; c < ARRAY_SIZE(cases); ++c) {
        const extreme_t* e = &cases[c];
        char* src = NULL;
        md_array_push_array(src, e->head, strlen(e->head), scratch);
        for (int i = 0; i < e->count; ++i) md_array_push_array(src, e->unit, strlen(e->unit), scratch);
        md_array_push_array(src, e->tail, strlen(e->tail), scratch);
        if (e->close) for (int i = 0; i < e->count; ++i) md_array_push(src, e->close, scratch);
        md_array_push(src, '\n', scratch);
        const size_t len = md_array_size(src);

        const char* err = fuzz_check_copy(src, len, scratch);
        if (err) printf("%s: %s\n", e->name, err);
        EXPECT_TRUE(err == NULL);

        md_vial_program_t* prog = md_vial_compile((str_t){src, len}, NULL, scratch);
        const md_vial_diagnostic_t* d = NULL;
        const size_t n = md_vial_diagnostics(prog, &d);
        size_t deep = 0;
        for (size_t i = 0; i < n; ++i) deep += strcmp(d[i].code, "E0104") == 0;
        if (e->deep && deep != 1) printf("%s: %zu E0104, not 1\n", e->name, deep);
        if (e->deep) EXPECT_EQ(1, deep);
        else         EXPECT_EQ(0, deep);
        md_vial_program_free(prog);
        md_arena_allocator_reset(scratch);
    }
    md_arena_allocator_destroy(scratch);
}

UTEST_F(vial_fuzz, mutations) {
    const size_t iterations = (size_t)fuzz_env("VIAL_FUZZ_ITERATIONS", 3000);
    fuzz_rng_t rng = {fuzz_env("VIAL_FUZZ_SEED", 0x5eed5eedull) | 1};
    md_allocator_i* scratch = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    const str_t* seeds = utest_fixture->seeds;
    const size_t num_seeds = md_array_size(seeds);

    size_t num_bytes = 0, num_failed = 0;
    for (size_t it = 0; it < iterations && num_failed < 3; ++it) {
        const str_t s = seeds[fuzz_below(&rng, num_seeds)];
        fuzz_buf_t b = {md_alloc(scratch, s.len + 1), s.len};
        memcpy(b.ptr, s.ptr, s.len);
        fuzz_mutate(&b, seeds, utest_fixture->dict, &rng, scratch);
        num_bytes += b.len;

        const char* err = fuzz_check_copy(b.ptr, b.len, scratch);
        if (err) {
            num_failed += 1;
            printf("vial_fuzz iteration %zu: %s\n", it, err);
            fuzz_print_input(b.ptr, b.len);
        }
        md_arena_allocator_reset(scratch);
    }
    printf("vial_fuzz: %zu mutated inputs, %zu bytes, %zu failed\n", iterations, num_bytes, num_failed);
    EXPECT_EQ(0, num_failed);
    md_arena_allocator_destroy(scratch);
}

#include "utest.h"

#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#include <vial/vial_internal.h>

#include <stdio.h>
#include <string.h>

// ### LEXER ###

typedef struct lexed_t {
    vial_token_t* tokens;
    vial_diags_t  diags;
    str_t         dump;
} lexed_t;

static lexed_t lex(const char* src, md_allocator_i* alloc) {
    lexed_t r = {0};
    const str_t s = str_from_cstr(src);
    vial_diags_init(&r.diags, s, alloc);
    r.tokens = vial_lex(s, &r.diags, alloc);
    r.dump = vial_tokens_dump(r.tokens, s, alloc);
    return r;
}

// The token kinds in order, space separated, with operators and keywords as their text
static str_t kinds(const lexed_t* r, const char* src, md_allocator_i* alloc) {
    char buf[2048] = {0};
    size_t len = 0;
    for (size_t i = 0; i < md_array_size(r->tokens); ++i) {
        const vial_token_t* t = &r->tokens[i];
        const char* k = vial_token_kind_str(t->kind);
        char item[128];
        if (t->kind == MD_VIAL_TOKEN_OPERATOR || t->kind == MD_VIAL_TOKEN_PUNCT) {
            snprintf(item, sizeof(item), "%s", vial_op_str((vial_op_t)t->sub));
        } else if (t->kind == MD_VIAL_TOKEN_KEYWORD) {
            snprintf(item, sizeof(item), "%s", vial_kw_str((vial_kw_t)t->sub));
        } else if (t->kind == MD_VIAL_TOKEN_NAME || t->kind == MD_VIAL_TOKEN_INT || t->kind == MD_VIAL_TOKEN_FLOAT) {
            snprintf(item, sizeof(item), "%s:%.*s", k, (int)(t->span.end - t->span.beg), src + t->span.beg);
        } else {
            snprintf(item, sizeof(item), "%s", k);
        }
        len += snprintf(buf + len, sizeof(buf) - len, "%s%s", len ? " " : "", item);
    }
    return str_copy_cstrn(buf, len, alloc);
}

#define EXPECT_KINDS(src, expected) do { \
    lexed_t r_ = lex(src, alloc); \
    str_t k_ = kinds(&r_, src, alloc); \
    if (!str_eq_cstr(k_, expected)) printf("got:      %.*s\nexpected: %s\n", (int)k_.len, k_.ptr, expected); \
    EXPECT_TRUE(str_eq_cstr(k_, expected)); \
    EXPECT_EQ(0, md_array_size(r_.diags.items)); \
} while (0)

#define EXPECT_DIAG(src, want_) do { \
    lexed_t r_ = lex(src, alloc); \
    const size_t n_ = vial_diag_count(&r_.diags, want_); \
    if (n_ != 1) { \
        printf("%s: expected one %s, got %zu of %zu\n", src, want_, n_, md_array_size(r_.diags.items)); \
        for (size_t i_ = 0; i_ < md_array_size(r_.diags.items); ++i_) \
            printf("  %s %.*s\n", r_.diags.items[i_].code, (int)r_.diags.items[i_].message.len, r_.diags.items[i_].message.ptr); \
    } \
    EXPECT_EQ(1, n_); \
} while (0)

struct vial_lex {
    md_allocator_i* alloc;
};

UTEST_F_SETUP(vial_lex) {
    utest_fixture->alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
}

UTEST_F_TEARDOWN(vial_lex) {
    md_arena_allocator_destroy(utest_fixture->alloc);
}

UTEST_F(vial_lex, names_and_keywords) {
    md_allocator_i* alloc = utest_fixture->alloc;
    EXPECT_KINDS("a = b", "NAME:a = NAME:b NEWLINE EOF");
    EXPECT_KINDS("x is None and not y", "NAME:x is None and not NAME:y NEWLINE EOF");
    EXPECT_KINDS("_core2 = residues", "NAME:_core2 = NAME:residues NEWLINE EOF");
    EXPECT_KINDS("match = 1", "NAME:match = INT:1 NEWLINE EOF");   // soft keywords are names
    EXPECT_KINDS("dö = 1", "NAME:dö = INT:1 NEWLINE EOF");         // Python's names are Unicode
}

UTEST_F(vial_lex, numbers_and_ranges) {
    md_allocator_i* alloc = utest_fixture->alloc;
    EXPECT_KINDS("1_000_000", "INT:1_000_000 NEWLINE EOF");
    EXPECT_KINDS("0.5 1e-3 2. .5", "FLOAT:0.5 FLOAT:1e-3 FLOAT:2. FLOAT:.5 NEWLINE EOF");
    EXPECT_KINDS("0x1F 0b101", "INT:0x1F INT:0b101 NEWLINE EOF");
    // A number followed by `..` ends there: 1..4 is a range, never 1. and .4
    EXPECT_KINDS("1..4", "INT:1 .. INT:4 NEWLINE EOF");
    EXPECT_KINDS("0.0..5.0", "FLOAT:0.0 .. FLOAT:5.0 NEWLINE EOF");
    EXPECT_KINDS("1..n+1", "INT:1 .. NAME:n + INT:1 NEWLINE EOF");
    EXPECT_KINDS("c[..., 0]", "NAME:c [ ... , INT:0 ] NEWLINE EOF");
    EXPECT_KINDS("a.com.x", "NAME:a . NAME:com . NAME:x NEWLINE EOF");
}

UTEST_F(vial_lex, operators) {
    md_allocator_i* alloc = utest_fixture->alloc;
    EXPECT_KINDS("a @ b ** c // d", "NAME:a @ NAME:b ** NAME:c // NAME:d NEWLINE EOF");
    EXPECT_KINDS("x += 1", "NAME:x += INT:1 NEWLINE EOF");
    EXPECT_KINDS("def f(a: group) -> float: pass",
        "def NAME:f ( NAME:a : NAME:group ) -> NAME:float : pass NEWLINE EOF");
    EXPECT_KINDS("2 < n <= 500 != m", "INT:2 < NAME:n <= INT:500 != NAME:m NEWLINE EOF");
}

UTEST_F(vial_lex, strings) {
    md_allocator_i* alloc = utest_fixture->alloc;
    EXPECT_KINDS("'CA' \"GTP\"", "STRING STRING NEWLINE EOF");
    EXPECT_KINDS("\"\"\"Doc\nstring\"\"\"", "STRING NEWLINE EOF");
    EXPECT_KINDS("r'\\d' 'it\\'s'", "STRING STRING NEWLINE EOF");
    EXPECT_KINDS("x = '#not a comment'", "NAME:x = STRING NEWLINE EOF");

    lexed_t r = lex("f'{x}'", alloc);
    EXPECT_EQ(MD_VIAL_TOKEN_STRING, r.tokens[0].kind);
    EXPECT_TRUE(r.tokens[0].flags & VIAL_STR_FORMAT);   // the parser refuses it
}

UTEST_F(vial_lex, lines_and_indentation) {
    md_allocator_i* alloc = utest_fixture->alloc;
    EXPECT_KINDS("def f(x):\n    return x\ny = 1\n",
        "def NAME:f ( NAME:x ) : NEWLINE INDENT return NAME:x NEWLINE DEDENT NAME:y = INT:1 NEWLINE EOF");
    // Blank and comment lines are not logical lines, whatever their indentation
    EXPECT_KINDS("if a:\n    b\n\n  # note\n    c\n",
        "if NAME:a : NEWLINE INDENT NAME:b NEWLINE COMMENT NAME:c NEWLINE DEDENT EOF");
    // Line breaks inside brackets, and after a backslash, join lines
    EXPECT_KINDS("x = f(1,\n      2)\n", "NAME:x = NAME:f ( INT:1 , INT:2 ) NEWLINE EOF");
    EXPECT_KINDS("x = 1 + \\\n    2\n", "NAME:x = INT:1 + INT:2 NEWLINE EOF");
    // Nested blocks close at the end of the text
    EXPECT_KINDS("for r in s:\n  if r:\n    pass",
        "for NAME:r in NAME:s : NEWLINE INDENT if NAME:r : NEWLINE INDENT pass NEWLINE DEDENT DEDENT EOF");
    // \r\n and a byte order mark
    EXPECT_KINDS("\xEF\xBB\xBF" "a = 1\r\nb = 2\r\n", "NAME:a = INT:1 NEWLINE NAME:b = INT:2 NEWLINE EOF");
    EXPECT_KINDS("x = 1  # -> int", "NAME:x = INT:1 COMMENT NEWLINE EOF");
}

UTEST_F(vial_lex, errors) {
    md_allocator_i* alloc = utest_fixture->alloc;
    EXPECT_DIAG("x = $", "E0101");
    EXPECT_DIAG("x = )", "E0101");
    EXPECT_DIAG("x = (1, 2]", "E0101");
    EXPECT_DIAG("x = 1abc", "E0101");
    EXPECT_DIAG("x = 1__0", "E0101");
    EXPECT_DIAG("x = 'abc\ny = 2", "E0102");
    EXPECT_DIAG("x = \"\"\"abc", "E0102");
    EXPECT_DIAG("x = f(1, 2\n", "E0102");
    EXPECT_DIAG("if a:\n    b\n   c\n", "E0103");
    EXPECT_DIAG("x = 1j", "E0803");

    // A message names the line the bracket was opened on
    lexed_t r = lex("a = 1\nx = f(1,\n", alloc);
    ASSERT_EQ(1, md_array_size(r.diags.items));
    EXPECT_TRUE(str_eq_cstr(r.diags.items[0].message, "`(` opened on line 2 is never closed"));

    // E0103 says how far the line is indented and what the block uses
    r = lex("if a:\n    b\n   c\n", alloc);
    ASSERT_EQ(1, md_array_size(r.diags.items));
    EXPECT_TRUE(str_eq_cstr(r.diags.items[0].message, "this line is indented by 3 spaces; the block uses 4"));
}

// A missing `)` does not swallow the statements after it
UTEST_F(vial_lex, recovery) {
    md_allocator_i* alloc = utest_fixture->alloc;
    lexed_t r = lex("a = f(1, 2\nb = 3\n", alloc);
    EXPECT_EQ(1, vial_diag_count(&r.diags, "E0102"));
    str_t k = kinds(&r, "a = f(1, 2\nb = 3\n", alloc);
    EXPECT_TRUE(str_eq_cstr(k, "NAME:a = NAME:f ( INT:1 , INT:2 NEWLINE NAME:b = INT:3 NEWLINE EOF"));

    // but a continued argument list on the next line stays joined
    EXPECT_KINDS("a = f(1,\n    b=2)\n", "NAME:a = NAME:f ( INT:1 , NAME:b = INT:2 ) NEWLINE EOF");
}

// Every byte of the text is covered by a token or is whitespace, and spans never go backwards
UTEST_F(vial_lex, spans) {
    md_allocator_i* alloc = utest_fixture->alloc;
    const char* src = "def f(a: group, b=None):\n    \"\"\"Doc.\"\"\"\n    return a.com - b  # c\n";
    lexed_t r = lex(src, alloc);
    EXPECT_EQ(0, md_array_size(r.diags.items));
    uint32_t prev = 0;
    for (size_t i = 0; i < md_array_size(r.tokens); ++i) {
        const md_vial_span_t s = r.tokens[i].span;
        EXPECT_LE(prev, s.beg);
        EXPECT_LE(s.beg, s.end);
        EXPECT_LE(s.end, (uint32_t)strlen(src));
        for (uint32_t k = prev; k < s.beg; ++k) {
            const char c = src[k];
            EXPECT_TRUE(c == ' ' || c == '\n' || c == '\t');
        }
        if (s.end > prev) prev = s.end;
    }
}

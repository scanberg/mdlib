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

// ### PARSER ###

#include <vial/ast.h>

typedef struct parsed_t {
    vial_tree_t  tree;
    vial_diags_t diags;
    str_t        src;
} parsed_t;

static parsed_t parse(const char* src, md_allocator_i* alloc) {
    parsed_t r = {.src = str_from_cstr(src)};
    vial_diags_init(&r.diags, r.src, alloc);
    vial_token_t* tokens = vial_lex(r.src, &r.diags, alloc);
    r.tree = vial_parse(tokens, r.src, &r.diags, alloc);
    return r;
}

// The tree as an s-expression without spans: (assign (name x) (int 1))
static void sexpr(char* buf, size_t cap, size_t* len, const parsed_t* r, vial_node_id id) {
#define OUT(...) do { int n_ = snprintf(buf + *len, cap - *len, __VA_ARGS__); if (n_ > 0) *len += (size_t)n_; if (*len >= cap) *len = cap - 1; } while (0)
    if (id == 0) { OUT("-"); return; }
    const vial_node_t* n = &r->tree.nodes[id];
    OUT("(%s", vial_node_kind_str((vial_node_kind_t)n->kind));
    switch (n->kind) {
    case VIAL_N_BINARY: case VIAL_N_UNARY: case VIAL_N_AUG_ASSIGN: OUT(" %s", vial_op_str((vial_op_t)n->sub)); break;
    case VIAL_N_BOOLOP: case VIAL_N_CONST: OUT(" %s", vial_kw_str((vial_kw_t)n->sub)); break;
    case VIAL_N_CMPOP: OUT(" %s", vial_cmp_str((vial_cmp_t)n->sub)); break;
    case VIAL_N_NAME: case VIAL_N_INT: case VIAL_N_FLOAT: case VIAL_N_STRING:
        OUT(" %.*s", (int)(n->span.end - n->span.beg), r->src.ptr + n->span.beg); break;
    default: break;
    }
    for (uint32_t i = 0; i < n->child_count; ++i) {
        OUT(" ");
        sexpr(buf, cap, len, r, r->tree.children[n->child_beg + i]);
    }
    OUT(")");
#undef OUT
}

static str_t tree_str(const parsed_t* r, md_allocator_i* alloc) {
    char buf[4096];
    size_t len = 0;
    buf[0] = 0;
    // The module's statements, without the module itself
    const vial_node_t* m = &r->tree.nodes[r->tree.root];
    for (uint32_t i = 0; i < m->child_count; ++i) {
        if (i) len += (size_t)snprintf(buf + len, sizeof(buf) - len, " ");
        sexpr(buf, sizeof(buf), &len, r, r->tree.children[m->child_beg + i]);
    }
    return str_copy_cstrn(buf, len, alloc);
}

#define EXPECT_TREE(src, expected) do { \
    parsed_t r_ = parse(src, alloc); \
    str_t t_ = tree_str(&r_, alloc); \
    if (!str_eq_cstr(t_, expected)) printf("source:   %s\ngot:      %.*s\nexpected: %s\n", src, (int)t_.len, t_.ptr, expected); \
    EXPECT_TRUE(str_eq_cstr(t_, expected)); \
    for (size_t i_ = 0; i_ < md_array_size(r_.diags.items); ++i_) \
        printf("unexpected %s %.*s\n", r_.diags.items[i_].code, (int)r_.diags.items[i_].message.len, r_.diags.items[i_].message.ptr); \
    EXPECT_EQ(0, md_array_size(r_.diags.items)); \
} while (0)

// Exactly these codes, in this order
#define EXPECT_CODES(src, expected) do { \
    parsed_t r_ = parse(src, alloc); \
    char got_[256] = {0}; size_t gl_ = 0; \
    for (size_t i_ = 0; i_ < md_array_size(r_.diags.items); ++i_) \
        gl_ += (size_t)snprintf(got_ + gl_, sizeof(got_) - gl_, "%s%s", i_ ? " " : "", r_.diags.items[i_].code); \
    if (strcmp(got_, expected) != 0) { \
        printf("source:   %s\ngot:      %s\nexpected: %s\n", src, got_, expected); \
        for (size_t i_ = 0; i_ < md_array_size(r_.diags.items); ++i_) \
            printf("  %s %.*s\n", r_.diags.items[i_].code, (int)r_.diags.items[i_].message.len, r_.diags.items[i_].message.ptr); \
    } \
    EXPECT_STREQ(expected, got_); \
} while (0)

struct vial_parse {
    md_allocator_i* alloc;
};

UTEST_F_SETUP(vial_parse) {
    utest_fixture->alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
}

UTEST_F_TEARDOWN(vial_parse) {
    md_arena_allocator_destroy(utest_fixture->alloc);
}

UTEST_F(vial_parse, precedence) {
    md_allocator_i* alloc = utest_fixture->alloc;
    EXPECT_TREE("a + b * c", "(expr (binary + (name a) (binary * (name b) (name c))))");
    EXPECT_TREE("-a ** 2", "(expr (unary - (binary ** (name a) (int 2))))");
    EXPECT_TREE("a ** -b", "(expr (binary ** (name a) (unary - (name b))))");
    EXPECT_TREE("a @ b + c", "(expr (binary + (binary @ (name a) (name b)) (name c)))");
    EXPECT_TREE("1..n+1", "(expr (range (int 1) (binary + (name n) (int 1))))");
    EXPECT_TREE("x in 10..20", "(expr (compare (name x) (cmpop in) (range (int 10) (int 20))))");
    EXPECT_TREE("2 < n <= 500", "(expr (compare (int 2) (cmpop <) (name n) (cmpop <=) (int 500)))");
    EXPECT_TREE("a is not None", "(expr (compare (name a) (cmpop is not) (const None)))");
    EXPECT_TREE("x not in s", "(expr (compare (name x) (cmpop not in) (name s)))");
    EXPECT_TREE("not a and b or c", "(expr (boolop or (boolop and (not (name a)) (name b)) (name c)))");
    EXPECT_TREE("a if c else b", "(expr (ternary (name a) (name c) (name b)))");
    EXPECT_TREE("a & b | c ^ d", "(expr (binary | (binary & (name a) (name b)) (binary ^ (name c) (name d))))");
    EXPECT_TREE("a and b and c", "(expr (boolop and (name a) (name b) (name c)))");
}

UTEST_F(vial_parse, primaries) {
    md_allocator_i* alloc = utest_fixture->alloc;
    EXPECT_TREE("f(a, k=1)", "(expr (call (name f) (name a) (kwarg (name k) (int 1))))");
    EXPECT_TREE("res_a.com.x", "(expr (attr (attr (name res_a) (name com)) (name x)))");
    EXPECT_TREE("c[..., 0]", "(expr (subscript (name c) (tuple (ellipsis) (int 0))))");
    EXPECT_TREE("c[:, None]", "(expr (subscript (name c) (tuple (slice - - -) (const None))))");
    EXPECT_TREE("a[1:]", "(expr (subscript (name a) (slice (int 1) - -)))");
    EXPECT_TREE("a[::2]", "(expr (subscript (name a) (slice - - (int 2))))");
    EXPECT_TREE("a[-3:]", "(expr (subscript (name a) (slice (unary - (int 3)) - -)))");
    EXPECT_TREE("v[:, [-1, -2]]", "(expr (subscript (name v) (tuple (slice - - -) (list (unary - (int 1)) (unary - (int 2))))))");
    EXPECT_TREE("x = ()", "(assign (name x) (tuple))");
    EXPECT_TREE("x = (1,)", "(assign (name x) (tuple (int 1)))");
    EXPECT_TREE("x = (1)", "(assign (name x) (int 1))");
    EXPECT_TREE("s = 'a' \"b\"", "(assign (name s) (string 'a' \"b\"))");
    EXPECT_TREE("[r for r in res if r]", "(expr (listcomp (name r) (compfor (name r) (name res) (name r))))");
    EXPECT_TREE("[x for x, y in zip(a, b)]",
        "(expr (listcomp (name x) (compfor (tuple (name x) (name y)) (call (name zip) (name a) (name b)))))");
    EXPECT_TREE("[[a for a in r] for r in s]",
        "(expr (listcomp (listcomp (name a) (compfor (name a) (name r))) (compfor (name r) (name s))))");
    EXPECT_TREE("[x for _ in range(n + 1)]",
        "(expr (listcomp (name x) (compfor (name _) (call (name range) (binary + (name n) (int 1))))))");
}

UTEST_F(vial_parse, statements) {
    md_allocator_i* alloc = utest_fixture->alloc;
    EXPECT_TREE("a, b = common(x, y)", "(assign (tuple (name a) (name b)) (call (name common) (name x) (name y)))");
    EXPECT_TREE("a = b = 1", "(assign (name a) (name b) (int 1))");
    EXPECT_TREE("x += 1", "(augassign += (name x) (int 1))");
    EXPECT_TREE("a = 1; b = 2", "(assign (name a) (int 1)) (assign (name b) (int 2))");
    EXPECT_TREE("def f(a: group, b, op='x') -> float:\n    return a\n",
        "(def (name f) (params (param (name a) (name group) -) (param (name b) - -) (param (name op) - (string 'x'))) "
        "(name float) (block (return (name a))))");
    EXPECT_TREE("if a:\n    pass\nelif b:\n    pass\nelse:\n    x = 1\n",
        "(if (name a) (block (pass)) (if (name b) (block (pass)) (block (assign (name x) (int 1)))))");
    EXPECT_TREE("for r in res:\n    if r: break\n",
        "(for (name r) (name res) (block (if (name r) (block (break)) -)))");
    EXPECT_TREE("assert op in ('a', 'b'), \"msg\"",
        "(assert (compare (name op) (cmpop in) (tuple (string 'a') (string 'b'))) (string \"msg\"))");
    EXPECT_TREE("def f():\n    def g(s, v):\n        return s + v\n    return g\n",
        "(def (name f) (params) - (block (def (name g) (params (param (name s) - -) (param (name v) - -)) - "
        "(block (return (binary + (name s) (name v))))) (return (name g))))");
    EXPECT_TREE("def f(x):\n    \"\"\"Doc.\"\"\"\n    return x\n",
        "(def (name f) (params (param (name x) - -)) - (block (expr (string \"\"\"Doc.\"\"\")) (return (name x))))");
}

UTEST_F(vial_parse, refused) {
    md_allocator_i* alloc = utest_fixture->alloc;
    EXPECT_CODES("class Foo:\n    pass\n", "E0803");
    EXPECT_CODES("import numpy", "E0803");
    EXPECT_CODES("from x import y", "E0803");
    EXPECT_CODES("f = lambda s, v: s + v", "E0803");
    EXPECT_CODES("try:\n    a = 1\nexcept:\n    pass\nfinally:\n    pass\n", "E0803");
    EXPECT_CODES("x = {1, 5, 9}", "E0803");
    EXPECT_CODES("x = {'a': 1}", "E0803");
    EXPECT_CODES("x = {}", "E0803");
    EXPECT_CODES("def f(*args): pass", "E0803");
    EXPECT_CODES("f(**k)", "E0803");
    EXPECT_CODES("x = f'{a}'", "E0803");
    EXPECT_CODES("x = (a for a in b)", "E0803");
    EXPECT_CODES("x = sum(a for a in b)", "E0803");
    EXPECT_CODES("global x", "E0803");
    EXPECT_CODES("@dec\ndef f(): pass", "E0803");
    EXPECT_CODES("with a as b:\n    pass\n", "E0803");
    EXPECT_CODES("del x", "E0803");
    EXPECT_CODES("x: int = 1", "E0803");
    EXPECT_CODES("for a in b:\n    pass\nelse:\n    pass\n", "E0803");
    EXPECT_CODES("if (n := 1): pass", "E0803");
    EXPECT_CODES("raise x", "E0803");
    EXPECT_CODES("def f():\n    yield 1\n", "E0803");
    // Not supported yet
    EXPECT_CODES("while x > 1.0:\n    x = x / 2.0\n", "E0801");
    EXPECT_CODES("def f(a):\n    a[0] = 1\n", "E0801");
    // v1 spellings
    EXPECT_CODES("residue(1:4)", "E0701");
    EXPECT_CODES("d = distance(1, 2) in residue(:)", "E0701");
    // A refused construct does not hide the next statement's error
    EXPECT_CODES("import x\ny = = 1\n", "E0803 E0101");
}

UTEST_F(vial_parse, errors_and_recovery) {
    md_allocator_i* alloc = utest_fixture->alloc;
    EXPECT_CODES("x = = 1", "E0101");
    EXPECT_CODES("x = 1 2", "E0101");
    EXPECT_CODES("if a\n    pass\n", "E0101");
    EXPECT_CODES("if a:\npass\n", "E0103");
    EXPECT_CODES("x = 1\n    y = 2\n", "E0103");
    EXPECT_CODES("def (a): pass", "E0101");
    EXPECT_CODES("x = a.1", "E0101");
    EXPECT_CODES("x = (1 + )", "E0101");
    // One mistake, one diagnostic: an unclosed bracket is the lexer's, not the parser's too
    EXPECT_CODES("a = f(1, 2\nb = 3\n", "E0102");
    EXPECT_CODES("a = (1, 2]\nb = 3\n", "E0101");
    EXPECT_CODES("x = $\ny = 1\n", "E0101");
    // Each broken statement is reported, and the ones between parse
    EXPECT_CODES("a = = 1\nb = 2\nc = * \nd = 4\n", "E0101 E0101");
    {
        parsed_t r = parse("a = = 1\nb = 2\nc = )\nd = 4\n", alloc);
        const vial_node_t* m = &r.tree.nodes[r.tree.root];
        EXPECT_EQ(4, m->child_count);   // every statement has a node, broken or not
    }
    // An error in a header still parses its block
    EXPECT_CODES("if a b:\n    x = = 1\n", "E0101 E0101");
}

// The tree's spans nest: every child lies within its parent
UTEST_F(vial_parse, spans_nest) {
    md_allocator_i* alloc = utest_fixture->alloc;
    const char* src = "def f(a: group, b=None):\n    d = norm(displacement(positions(a)[:, None], positions(b)))\n"
                      "    return min(d) if b is None else max(d)\nx = [r.com for r in residues(1..4) if r]\n";
    parsed_t r = parse(src, alloc);
    EXPECT_EQ(0, md_array_size(r.diags.items));
    for (size_t id = 1; id < md_array_size(r.tree.nodes); ++id) {
        const vial_node_t* n = &r.tree.nodes[id];
        EXPECT_LE(n->span.beg, n->span.end);
        for (uint32_t i = 0; i < n->child_count; ++i) {
            const vial_node_id c = r.tree.children[n->child_beg + i];
            if (!c) continue;
            EXPECT_LE(n->span.beg, r.tree.nodes[c].span.beg);
            EXPECT_LE(r.tree.nodes[c].span.end, n->span.end);
        }
    }
}

// ### GOLDEN DUMPS ###
// unittest/vial_golden/<name>.vial is parsed and its tree and diagnostics compared with <name>.out. A change to
// either is a reviewed diff. Run with VIAL_UPDATE_GOLDEN=1 to write the .out files anew.

#include <stdlib.h>

#define VIAL_GOLDEN_DIR MD_UNITTEST_SOURCE_DIR "/vial_golden/"

static str_t golden_output(str_t src, md_allocator_i* alloc) {
    vial_diags_t diags;
    vial_diags_init(&diags, src, alloc);
    vial_token_t* tokens = vial_lex(src, &diags, alloc);
    vial_tree_t tree = vial_parse(tokens, src, &diags, alloc);
    str_t dump = vial_tree_dump(&tree, src, alloc);

    char* buf = NULL;
    md_array_push_array(buf, dump.ptr, dump.len, alloc);
    const char* sep = "--- diagnostics\n";
    md_array_push_array(buf, sep, strlen(sep), alloc);
    for (size_t i = 0; i < md_array_size(diags.items); ++i) {
        const md_vial_diagnostic_t* d = &diags.items[i];
        const uint32_t line = vial_line_of(&diags, d->span.beg);
        const uint32_t col = d->span.beg - diags.line_beg[line - 1] + 1;
        char tmp[640];
        int n = snprintf(tmp, sizeof(tmp), "%s %u:%u %.*s%s%.*s\n", d->code, line, col, (int)d->message.len, d->message.ptr,
            d->fix.len ? " | fix: " : "", (int)d->fix.len, d->fix.ptr);
        if (n > 0) md_array_push_array(buf, tmp, (size_t)(n < (int)sizeof(tmp) ? n : (int)sizeof(tmp) - 1), alloc);
    }
    return (str_t){buf, md_array_size(buf)};
}

UTEST_F(vial_parse, golden) {
    md_allocator_i* alloc = utest_fixture->alloc;
    static const char* names[] = {"expressions", "statements", "refused", "errors"};
    const bool update = getenv("VIAL_UPDATE_GOLDEN") != NULL;
    for (size_t i = 0; i < ARRAY_SIZE(names); ++i) {
        char in_path[512], out_path[512];
        snprintf(in_path, sizeof(in_path), VIAL_GOLDEN_DIR "%s.vial", names[i]);
        snprintf(out_path, sizeof(out_path), VIAL_GOLDEN_DIR "%s.out", names[i]);
        str_t src = load_textfile(str_from_cstr(in_path), alloc);
        ASSERT_FALSE(str_empty(src));
        str_t got = golden_output(src, alloc);
        if (update) {
            FILE* f = fopen(out_path, "wb");
            ASSERT_TRUE(f != NULL);
            fwrite(got.ptr, 1, got.len, f);
            fclose(f);
            continue;
        }
        str_t want = load_textfile(str_from_cstr(out_path), alloc);
        const bool same = str_eq(got, want);
        if (!same) printf("%s differs from its golden dump; run with VIAL_UPDATE_GOLDEN=1 and review the diff\n", out_path);
        EXPECT_TRUE(same);
    }
}

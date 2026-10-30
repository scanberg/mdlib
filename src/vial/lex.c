// The VIAL lexer: Python 3's tokens, plus the inclusive range `..`.
//
// A logical line ends in NEWLINE; a line holding only a comment or nothing is not a logical line. Inside brackets a
// line break joins lines, as in Python, and so does a backslash at the end of a line. Indentation becomes INDENT and
// DEDENT tokens. Every token has a byte span, comments included, so the editor can colour all of the text.
//
// Errors are reported and lexing goes on, so broken text still gives tokens:
//   E0101  a character or number that is not VIAL
//   E0102  a bracket or string that is never closed
//   E0103  indentation that matches no enclosing block, or mixes tabs and spaces ambiguously

#include "vial_internal.h"

#include <core/md_common.h>

#include <string.h>
#include <stdio.h>
#include <stdarg.h>

static const char* kw_names[VIAL_KW_COUNT] = {
    [VIAL_KW_FALSE] = "False", [VIAL_KW_NONE_] = "None", [VIAL_KW_TRUE] = "True", [VIAL_KW_AND] = "and",
    [VIAL_KW_AS] = "as", [VIAL_KW_ASSERT] = "assert", [VIAL_KW_ASYNC] = "async", [VIAL_KW_AWAIT] = "await",
    [VIAL_KW_BREAK] = "break", [VIAL_KW_CLASS] = "class", [VIAL_KW_CONTINUE] = "continue", [VIAL_KW_DEF] = "def",
    [VIAL_KW_DEL] = "del", [VIAL_KW_ELIF] = "elif", [VIAL_KW_ELSE] = "else", [VIAL_KW_EXCEPT] = "except",
    [VIAL_KW_FINALLY] = "finally", [VIAL_KW_FOR] = "for", [VIAL_KW_FROM] = "from", [VIAL_KW_GLOBAL] = "global",
    [VIAL_KW_IF] = "if", [VIAL_KW_IMPORT] = "import", [VIAL_KW_IN] = "in", [VIAL_KW_IS] = "is",
    [VIAL_KW_LAMBDA] = "lambda", [VIAL_KW_NONLOCAL] = "nonlocal", [VIAL_KW_NOT] = "not", [VIAL_KW_OR] = "or",
    [VIAL_KW_PASS] = "pass", [VIAL_KW_RAISE] = "raise", [VIAL_KW_RETURN] = "return", [VIAL_KW_TRY] = "try",
    [VIAL_KW_WHILE] = "while", [VIAL_KW_WITH] = "with", [VIAL_KW_YIELD] = "yield",
};

static const char* op_names[VIAL_OP_COUNT] = {
    [VIAL_OP_LPAREN] = "(", [VIAL_OP_RPAREN] = ")", [VIAL_OP_LBRACKET] = "[", [VIAL_OP_RBRACKET] = "]",
    [VIAL_OP_LBRACE] = "{", [VIAL_OP_RBRACE] = "}", [VIAL_OP_COMMA] = ",", [VIAL_OP_COLON] = ":",
    [VIAL_OP_SEMICOLON] = ";", [VIAL_OP_DOT] = ".", [VIAL_OP_RANGE] = "..", [VIAL_OP_ELLIPSIS] = "...",
    [VIAL_OP_ARROW] = "->", [VIAL_OP_ASSIGN] = "=", [VIAL_OP_WALRUS] = ":=",
    [VIAL_OP_PLUS] = "+", [VIAL_OP_MINUS] = "-", [VIAL_OP_STAR] = "*", [VIAL_OP_SLASH] = "/", [VIAL_OP_DSLASH] = "//",
    [VIAL_OP_PERCENT] = "%", [VIAL_OP_POW] = "**", [VIAL_OP_AT] = "@", [VIAL_OP_LSHIFT] = "<<", [VIAL_OP_RSHIFT] = ">>",
    [VIAL_OP_AMP] = "&", [VIAL_OP_PIPE] = "|", [VIAL_OP_CARET] = "^", [VIAL_OP_TILDE] = "~",
    [VIAL_OP_LT] = "<", [VIAL_OP_GT] = ">", [VIAL_OP_LE] = "<=", [VIAL_OP_GE] = ">=", [VIAL_OP_EQ] = "==",
    [VIAL_OP_NE] = "!=",
    [VIAL_OP_PLUS_EQ] = "+=", [VIAL_OP_MINUS_EQ] = "-=", [VIAL_OP_STAR_EQ] = "*=", [VIAL_OP_SLASH_EQ] = "/=",
    [VIAL_OP_DSLASH_EQ] = "//=", [VIAL_OP_PERCENT_EQ] = "%=", [VIAL_OP_POW_EQ] = "**=", [VIAL_OP_AT_EQ] = "@=",
    [VIAL_OP_AMP_EQ] = "&=", [VIAL_OP_PIPE_EQ] = "|=", [VIAL_OP_CARET_EQ] = "^=", [VIAL_OP_LSHIFT_EQ] = "<<=",
    [VIAL_OP_RSHIFT_EQ] = ">>=",
};

const char* vial_kw_str(vial_kw_t kw) { return (kw > 0 && kw < VIAL_KW_COUNT) ? kw_names[kw] : "?"; }
const char* vial_op_str(vial_op_t op) { return (op > 0 && op < VIAL_OP_COUNT) ? op_names[op] : "?"; }

const char* vial_token_kind_str(md_vial_token_kind_t kind) {
    switch (kind) {
    case MD_VIAL_TOKEN_NAME:     return "NAME";
    case MD_VIAL_TOKEN_KEYWORD:  return "KW";
    case MD_VIAL_TOKEN_INT:      return "INT";
    case MD_VIAL_TOKEN_FLOAT:    return "FLOAT";
    case MD_VIAL_TOKEN_UNIT:     return "UNIT";
    case MD_VIAL_TOKEN_STRING:   return "STRING";
    case MD_VIAL_TOKEN_OPERATOR: return "OP";
    case MD_VIAL_TOKEN_PUNCT:    return "PUNCT";
    case MD_VIAL_TOKEN_COMMENT:  return "COMMENT";
    case MD_VIAL_TOKEN_NEWLINE:  return "NEWLINE";
    case MD_VIAL_TOKEN_INDENT:   return "INDENT";
    case MD_VIAL_TOKEN_DEDENT:   return "DEDENT";
    case MD_VIAL_TOKEN_ERROR:    return "ERROR";
    case MD_VIAL_TOKEN_EOF:      return "EOF";
    default: return "?";
    }
}

typedef struct indent_t {
    uint32_t col8;  // the column with tabs to the next multiple of 8
    uint32_t col1;  // the column with a tab counted as one, to tell an ambiguous mix of tabs and spaces
} indent_t;

typedef struct open_t {
    uint32_t pos;
    uint32_t line_col8;     // the indentation of the line it is on
    char     c;
} open_t;

typedef struct lexer_t {
    str_t           src;
    size_t          pos;
    vial_token_t*   tokens;
    open_t*         open;       // brackets that are open
    indent_t*       indents;    // the enclosing blocks' indentation; the first is 0
    bool            line_has_token;
    // The last answer of next_line_is_statement, kept for the blank lines that lead to the same line: one scan
    // over a run of blank lines inside brackets, not one per line
    size_t          ahead_from, ahead_to;
    uint32_t        ahead_col8;
    bool            ahead_result, ahead_valid;
    // Where the last opener's line starts, and its indentation, so that finding an opener's line costs time in
    // proportion to the text since the last opener, not to the length of its line
    size_t          col_scan_pos, col_line_beg;
    uint32_t        col_col8;
    bool            col_valid;
    vial_diags_t*   diags;
    md_allocator_i* alloc;
} lexer_t;

static inline int peek(const lexer_t* l, size_t ahead) {
    const size_t p = l->pos + ahead;
    return p < l->src.len ? (unsigned char)l->src.ptr[p] : 0;
}

static inline bool at_end(const lexer_t* l) { return l->pos >= l->src.len; }

static inline bool is_ident_start(int c) { return is_alpha(c) || c == '_' || c >= 0x80; }
static inline bool is_ident_char(int c)  { return is_ident_start(c) || is_digit(c); }

static inline md_vial_span_t span(size_t beg, size_t end) { return (md_vial_span_t){(uint32_t)beg, (uint32_t)end}; }

static void emit(lexer_t* l, md_vial_token_kind_t kind, size_t beg, size_t end, uint16_t sub, uint16_t flags) {
    vial_token_t tok = {.kind = kind, .span = span(beg, end), .sub = sub, .flags = flags};
    md_array_push(l->tokens, tok, l->alloc);
    if (kind != MD_VIAL_TOKEN_COMMENT && kind != MD_VIAL_TOKEN_NEWLINE && kind != MD_VIAL_TOKEN_INDENT &&
        kind != MD_VIAL_TOKEN_DEDENT) {
        l->line_has_token = true;
    }
}

static void emit_newline(lexer_t* l, size_t beg, size_t end) {
    if (l->line_has_token) emit(l, MD_VIAL_TOKEN_NEWLINE, beg, end, 0, 0);
    l->line_has_token = false;
}

static vial_kw_t keyword_of(str_t s) {
    for (int k = 1; k < VIAL_KW_COUNT; ++k) {
        if (str_eq_cstr(s, kw_names[k])) return (vial_kw_t)k;
    }
    return VIAL_KW_NONE;
}

static bool is_newline_at(const lexer_t* l, size_t p) {
    return p < l->src.len && (l->src.ptr[p] == '\n' || l->src.ptr[p] == '\r');
}

// The length of the line break at p: 2 for \r\n, 1 for \n or \r, 0 for none
static size_t newline_len(const lexer_t* l, size_t p) {
    if (p >= l->src.len) return 0;
    if (l->src.ptr[p] == '\r') return (p + 1 < l->src.len && l->src.ptr[p + 1] == '\n') ? 2 : 1;
    return l->src.ptr[p] == '\n' ? 1 : 0;
}

// Indentation of the line starting at p, as columns, and where its first non-blank byte is
static indent_t measure_indent(const lexer_t* l, size_t p, size_t* out_end) {
    indent_t ind = {0};
    while (p < l->src.len) {
        const char c = l->src.ptr[p];
        if (c == ' ')       { ind.col8 += 1; ind.col1 += 1; }
        else if (c == '\t') { ind.col8 = (ind.col8 / 8 + 1) * 8; ind.col1 += 1; }
        else if (c == '\f') { ind.col8 = 0; ind.col1 = 0; }
        else break;
        ++p;
    }
    *out_end = p;
    return ind;
}

// After a line break inside brackets: does the next line look like a new statement, at or left of the indentation
// of the line the outermost open bracket is on? Then the bracket was most likely never closed, and lexing goes on
// as if it had been, so one missing ')' does not swallow the rest of the script.
static bool next_line_is_statement_(const lexer_t* l, size_t p, uint32_t open_line_col8, size_t* line);

static bool next_line_is_statement(lexer_t* l, size_t p, uint32_t open_line_col8) {
    if (l->ahead_valid && l->ahead_col8 == open_line_col8 && l->ahead_from <= p && p <= l->ahead_to) return l->ahead_result;
    size_t line = l->src.len;
    const bool r = next_line_is_statement_(l, p, open_line_col8, &line);
    l->ahead_valid = true;
    l->ahead_from = p;
    l->ahead_to = line;
    l->ahead_col8 = open_line_col8;
    l->ahead_result = r;
    return r;
}

static bool next_line_is_statement_(const lexer_t* l, size_t p, uint32_t open_line_col8, size_t* line) {
    for (;;) {
        *line = p;
        size_t q;
        const indent_t ind = measure_indent(l, p, &q);
        if (q >= l->src.len) return false;
        const char c = l->src.ptr[q];
        if (c == '#' || c == '\n' || c == '\r') {   // blank or comment line: look at the next one
            while (q < l->src.len && !is_newline_at(l, q)) ++q;
            if (q >= l->src.len) return false;
            p = q + newline_len(l, q);
            continue;
        }
        if (ind.col8 > open_line_col8) return false;
        if (!is_ident_start((unsigned char)c)) return false;
        size_t e = q;
        while (e < l->src.len && is_ident_char((unsigned char)l->src.ptr[e])) ++e;
        const str_t word = {l->src.ptr + q, e - q};
        switch (keyword_of(word)) {
        case VIAL_KW_DEF: case VIAL_KW_CLASS: case VIAL_KW_IF: case VIAL_KW_ELIF: case VIAL_KW_ELSE:
        case VIAL_KW_FOR: case VIAL_KW_WHILE: case VIAL_KW_RETURN: case VIAL_KW_IMPORT: case VIAL_KW_FROM:
        case VIAL_KW_WITH: case VIAL_KW_TRY: case VIAL_KW_PASS: case VIAL_KW_BREAK: case VIAL_KW_CONTINUE:
        case VIAL_KW_ASSERT:
            return true;
        default: break;
        }
        // name = ..., but not name == ...
        while (e < l->src.len && (l->src.ptr[e] == ' ' || l->src.ptr[e] == '\t')) ++e;
        return e + 1 < l->src.len && l->src.ptr[e] == '=' && l->src.ptr[e + 1] != '=';
    }
}

// The indentation of the line an opener at pos is on. Openers come in order, so each call scans back only to
// the one before.
static uint32_t opener_line_col8(lexer_t* l, size_t pos) {
    const char* s = l->src.ptr;
    size_t p = pos;
    while (p > l->col_scan_pos && s[p - 1] != '\n' && s[p - 1] != '\r') --p;
    const bool same_line = l->col_valid && p == l->col_scan_pos && p > 0 && s[p - 1] != '\n' && s[p - 1] != '\r';
    const size_t line_beg = same_line ? l->col_line_beg : p;
    if (!l->col_valid || line_beg != l->col_line_beg) {
        size_t q;
        l->col_col8 = measure_indent(l, line_beg, &q).col8;
    }
    l->col_valid = true;
    l->col_scan_pos = pos;
    l->col_line_beg = line_beg;
    return l->col_col8;
}

static void close_all_open(lexer_t* l) {
    const open_t* o = &l->open[0];
    vial_diag(l->diags, "E0102", span(o->pos, o->pos + 1), "`%c` opened on line %u is never closed", o->c,
        vial_line_of(l->diags, o->pos));
    md_array_shrink(l->open, 0);
}

// The start of a logical line: compare its indentation with the enclosing blocks
static void handle_indent(lexer_t* l, indent_t ind, size_t p) {
    indent_t top = *md_array_last(l->indents);
    if (ind.col8 > top.col8) {
        if (ind.col1 <= top.col1) {
            vial_diag(l->diags, "E0103", span(p, p), "this line mixes tabs and spaces in a way that is ambiguous");
        }
        md_array_push(l->indents, ind, l->alloc);
        emit(l, MD_VIAL_TOKEN_INDENT, p, p, 0, 0);
        return;
    }
    if (ind.col8 == top.col8) {
        if (ind.col1 != top.col1) {
            vial_diag(l->diags, "E0103", span(p, p), "this line mixes tabs and spaces in a way that is ambiguous");
        }
        return;
    }
    // Close blocks down to the one this line belongs to. A line between two levels stays in the inner block,
    // which is what a typo of one space usually means, and is reported.
    while (md_array_size(l->indents) > 1 && ind.col8 < md_array_last(l->indents)->col8) {
        const indent_t outer = l->indents[md_array_size(l->indents) - 2];
        if (outer.col8 < ind.col8) {
            vial_diag(l->diags, "E0103", span(p, p), "this line is indented by %u spaces; the block uses %u",
                ind.col8, md_array_last(l->indents)->col8);
            return;
        }
        md_array_pop(l->indents);
        emit(l, MD_VIAL_TOKEN_DEDENT, p, p, 0, 0);
    }
}

static void lex_string(lexer_t* l, size_t beg, uint16_t flags) {
    const char q = l->src.ptr[l->pos];
    const bool triple = peek(l, 1) == q && peek(l, 2) == q;
    if (triple) flags |= VIAL_STR_TRIPLE;
    const size_t quote_pos = l->pos;
    l->pos += triple ? 3 : 1;
    for (;;) {
        if (at_end(l)) break;
        const char c = l->src.ptr[l->pos];
        if (c == '\\') {
            l->pos += 1;
            if (!at_end(l)) l->pos += newline_len(l, l->pos) ? newline_len(l, l->pos) : 1;
            continue;
        }
        if (!triple && (c == '\n' || c == '\r')) break;
        if (c == q) {
            if (!triple) { l->pos += 1; emit(l, MD_VIAL_TOKEN_STRING, beg, l->pos, 0, flags); return; }
            if (peek(l, 1) == q && peek(l, 2) == q) { l->pos += 3; emit(l, MD_VIAL_TOKEN_STRING, beg, l->pos, 0, flags); return; }
        }
        l->pos += 1;
    }
    const uint32_t line = vial_line_of(l->diags, (uint32_t)quote_pos);
    if (triple) {
        vial_diag(l->diags, "E0102", span(quote_pos, quote_pos + 3), "`%c%c%c` opened on line %u is never closed", q, q, q, line);
    } else {
        vial_diag(l->diags, "E0102", span(quote_pos, quote_pos + 1), "`%c` opened on line %u is never closed", q, line);
    }
    emit(l, MD_VIAL_TOKEN_STRING, beg, l->pos, 0, flags | VIAL_STR_UNCLOSED);
}

static size_t scan_digits(const lexer_t* l, size_t p, bool (*ok)(int), bool* bad_underscore) {
    // Digits with single underscores between them
    bool prev_digit = false;
    while (p < l->src.len) {
        const int c = (unsigned char)l->src.ptr[p];
        if (ok(c)) { prev_digit = true; ++p; continue; }
        if (c == '_') {
            const int n = p + 1 < l->src.len ? (unsigned char)l->src.ptr[p + 1] : 0;
            if (!prev_digit || !ok(n)) *bad_underscore = true;
            prev_digit = false;
            ++p;
            continue;
        }
        break;
    }
    return p;
}

static bool is_dec(int c) { return is_digit(c); }
static bool is_hex(int c) { return is_digit(c) || ('a' <= c && c <= 'f') || ('A' <= c && c <= 'F'); }
static bool is_oct(int c) { return '0' <= c && c <= '7'; }
static bool is_bin(int c) { return c == '0' || c == '1'; }

static void lex_number(lexer_t* l) {
    const size_t beg = l->pos;
    size_t p = l->pos;
    bool bad = false;
    bool is_float = false;
    const int c0 = peek(l, 0), c1 = peek(l, 1);

    if (c0 == '0' && (c1 == 'x' || c1 == 'X' || c1 == 'o' || c1 == 'O' || c1 == 'b' || c1 == 'B')) {
        bool (*ok)(int) = (c1 == 'x' || c1 == 'X') ? is_hex : (c1 == 'o' || c1 == 'O') ? is_oct : is_bin;
        const size_t d = p + 2;
        p = scan_digits(l, d, ok, &bad);
        if (p == d) bad = true;
    } else {
        p = scan_digits(l, p, is_dec, &bad);
        // A fraction, unless the dot starts `..`: 1..4 is a range, never 1. and .4
        if (p < l->src.len && l->src.ptr[p] == '.' && !(p + 1 < l->src.len && l->src.ptr[p + 1] == '.')) {
            p += 1;
            is_float = true;
            p = scan_digits(l, p, is_dec, &bad);
        }
        if (p < l->src.len && (l->src.ptr[p] == 'e' || l->src.ptr[p] == 'E')) {
            size_t e = p + 1;
            if (e < l->src.len && (l->src.ptr[e] == '+' || l->src.ptr[e] == '-')) ++e;
            if (e < l->src.len && is_digit((unsigned char)l->src.ptr[e])) {
                p = scan_digits(l, e, is_dec, &bad);
                is_float = true;
            }
        }
    }
    l->pos = p;

    if (!at_end(l) && (peek(l, 0) == 'j' || peek(l, 0) == 'J')) {
        l->pos += 1;
        vial_diag(l->diags, "E0803", span(beg, l->pos), "complex numbers are not part of VIAL");
        emit(l, MD_VIAL_TOKEN_ERROR, beg, l->pos, 0, 0);
        return;
    }
    if (!at_end(l) && is_ident_char(peek(l, 0))) {
        // A unit written against the number, 10ns, is lexed as the number and then the unit
        size_t e = l->pos;
        while (e < l->src.len && is_ident_char((unsigned char)l->src.ptr[e])) ++e;
        if (!bad && vial_unit_name((str_t){l->src.ptr + l->pos, e - l->pos}, NULL)) {
            emit(l, is_float ? MD_VIAL_TOKEN_FLOAT : MD_VIAL_TOKEN_INT, beg, l->pos, 0, 0);
            return;
        }
        l->pos = e;
        bad = true;
    }
    if (bad) {
        char ex[VIAL_EXCERPT_MAX + 4];
        const str_t e = vial_excerpt(ex, sizeof(ex), (str_t){l->src.ptr + beg, l->pos - beg});
        vial_diag(l->diags, "E0101", span(beg, l->pos), "`%.*s` is not a number", (int)e.len, e.ptr);
        emit(l, MD_VIAL_TOKEN_ERROR, beg, l->pos, 0, 0);
        return;
    }
    emit(l, is_float ? MD_VIAL_TOKEN_FLOAT : MD_VIAL_TOKEN_INT, beg, l->pos, 0, 0);
}

// Operators and punctuation, longest first
typedef struct op_entry_t { const char* str; uint8_t len; vial_op_t op; } op_entry_t;
static const op_entry_t op_table[] = {
    {"**=", 3, VIAL_OP_POW_EQ}, {"//=", 3, VIAL_OP_DSLASH_EQ}, {">>=", 3, VIAL_OP_RSHIFT_EQ},
    {"<<=", 3, VIAL_OP_LSHIFT_EQ}, {"...", 3, VIAL_OP_ELLIPSIS},
    {"**", 2, VIAL_OP_POW}, {"//", 2, VIAL_OP_DSLASH}, {"<<", 2, VIAL_OP_LSHIFT}, {">>", 2, VIAL_OP_RSHIFT},
    {"<=", 2, VIAL_OP_LE}, {">=", 2, VIAL_OP_GE}, {"==", 2, VIAL_OP_EQ}, {"!=", 2, VIAL_OP_NE},
    {"->", 2, VIAL_OP_ARROW}, {"+=", 2, VIAL_OP_PLUS_EQ}, {"-=", 2, VIAL_OP_MINUS_EQ}, {"*=", 2, VIAL_OP_STAR_EQ},
    {"/=", 2, VIAL_OP_SLASH_EQ}, {"%=", 2, VIAL_OP_PERCENT_EQ}, {"@=", 2, VIAL_OP_AT_EQ}, {"&=", 2, VIAL_OP_AMP_EQ},
    {"|=", 2, VIAL_OP_PIPE_EQ}, {"^=", 2, VIAL_OP_CARET_EQ}, {"..", 2, VIAL_OP_RANGE}, {":=", 2, VIAL_OP_WALRUS},
    {"+", 1, VIAL_OP_PLUS}, {"-", 1, VIAL_OP_MINUS}, {"*", 1, VIAL_OP_STAR}, {"/", 1, VIAL_OP_SLASH},
    {"%", 1, VIAL_OP_PERCENT}, {"@", 1, VIAL_OP_AT}, {"&", 1, VIAL_OP_AMP}, {"|", 1, VIAL_OP_PIPE},
    {"^", 1, VIAL_OP_CARET}, {"~", 1, VIAL_OP_TILDE}, {"<", 1, VIAL_OP_LT}, {">", 1, VIAL_OP_GT},
    {"(", 1, VIAL_OP_LPAREN}, {")", 1, VIAL_OP_RPAREN}, {"[", 1, VIAL_OP_LBRACKET}, {"]", 1, VIAL_OP_RBRACKET},
    {"{", 1, VIAL_OP_LBRACE}, {"}", 1, VIAL_OP_RBRACE}, {",", 1, VIAL_OP_COMMA}, {":", 1, VIAL_OP_COLON},
    {";", 1, VIAL_OP_SEMICOLON}, {".", 1, VIAL_OP_DOT}, {"=", 1, VIAL_OP_ASSIGN},
};

static bool is_punct(vial_op_t op) {
    switch (op) {
    case VIAL_OP_LPAREN: case VIAL_OP_RPAREN: case VIAL_OP_LBRACKET: case VIAL_OP_RBRACKET:
    case VIAL_OP_LBRACE: case VIAL_OP_RBRACE: case VIAL_OP_COMMA: case VIAL_OP_COLON: case VIAL_OP_SEMICOLON:
        return true;
    default: return false;
    }
}

static char closer_of(char open) { return open == '(' ? ')' : open == '[' ? ']' : '}'; }

static void lex_op(lexer_t* l) {
    const size_t beg = l->pos;
    for (size_t i = 0; i < ARRAY_SIZE(op_table); ++i) {
        const op_entry_t* e = &op_table[i];
        if (l->pos + e->len <= l->src.len && memcmp(l->src.ptr + l->pos, e->str, e->len) == 0) {
            const char c = e->str[0];
            if (e->len == 1 && (c == '(' || c == '[' || c == '{')) {
                open_t o = {.pos = (uint32_t)beg, .line_col8 = opener_line_col8(l, beg), .c = c};
                md_array_push(l->open, o, l->alloc);
            } else if (e->len == 1 && (c == ')' || c == ']' || c == '}')) {
                const size_t n = md_array_size(l->open);
                if (n == 0) {
                    l->pos += 1;
                    vial_diag(l->diags, "E0101", span(beg, l->pos), "unexpected `%c`", c);
                    emit(l, MD_VIAL_TOKEN_ERROR, beg, l->pos, e->op, 0);
                    return;
                }
                if (closer_of(l->open[n - 1].c) != c) {
                    // Close the nearest opener this matches, if any; otherwise the closer is stray
                    size_t k = n;
                    while (k > 0 && closer_of(l->open[k - 1].c) != c) --k;
                    const open_t inner = l->open[n - 1];
                    vial_diag(l->diags, "E0101", span(beg, beg + 1), "`%c` does not close `%c` opened on line %u",
                        c, inner.c, vial_line_of(l->diags, inner.pos));
                    // It closes the nearest opener it matches; when none does, it was most likely meant to close
                    // the innermost one. Either way it stands as a closer, marked as reported, so the parser takes
                    // it without a second diagnostic.
                    md_array_shrink(l->open, k > 0 ? k - 1 : n - 1);
                    l->pos += 1;
                    emit(l, MD_VIAL_TOKEN_PUNCT, beg, l->pos, (uint16_t)e->op, VIAL_TOK_REPORTED);
                    return;
                } else {
                    md_array_pop(l->open);
                }
            }
            l->pos += e->len;
            emit(l, is_punct(e->op) ? MD_VIAL_TOKEN_PUNCT : MD_VIAL_TOKEN_OPERATOR, beg, l->pos, (uint16_t)e->op, 0);
            return;
        }
    }
    // Not VIAL: one character, a whole UTF-8 sequence if it is one
    size_t len = 1;
    const unsigned char c = (unsigned char)l->src.ptr[l->pos];
    if (c >= 0xC0) while (l->pos + len < l->src.len && ((unsigned char)l->src.ptr[l->pos + len] & 0xC0) == 0x80) ++len;
    l->pos += len;
    vial_diag(l->diags, "E0101", span(beg, l->pos), "unexpected `%.*s`", (int)len, l->src.ptr + beg);
    emit(l, MD_VIAL_TOKEN_ERROR, beg, l->pos, 0, 0);
}

// A string prefix: r, b, u, f and the pairs rb, br, rf, fr, in either case
static bool string_prefix(str_t s, uint16_t* flags) {
    if (s.len == 0 || s.len > 2) return false;
    uint16_t f = 0;
    for (size_t i = 0; i < s.len; ++i) {
        const int c = to_lower(s.ptr[i]);
        uint16_t bit = c == 'r' ? VIAL_STR_RAW : c == 'b' ? VIAL_STR_BYTES : c == 'f' ? VIAL_STR_FORMAT : c == 'u' ? 0x100 : 0;
        if (!bit || (f & bit)) return false;
        f |= bit;
    }
    if ((f & 0x100) && s.len > 1) return false;
    if ((f & VIAL_STR_BYTES) && (f & VIAL_STR_FORMAT)) return false;
    *flags = f & (uint16_t)~0x100;
    return true;
}

vial_token_t* vial_lex(str_t src, vial_diags_t* diags, md_allocator_i* alloc) {
    lexer_t l = {.src = src, .diags = diags, .alloc = alloc};
    md_array_push(l.indents, ((indent_t){0, 0}), alloc);

    // A byte order mark is not text
    if (src.len >= 3 && (unsigned char)src.ptr[0] == 0xEF && (unsigned char)src.ptr[1] == 0xBB &&
        (unsigned char)src.ptr[2] == 0xBF) {
        l.pos = 3;
    }

    bool at_line_start = true;
    while (!at_end(&l)) {
        if (at_line_start && md_array_size(l.open) == 0) {
            size_t p;
            const indent_t ind = measure_indent(&l, l.pos, &p);
            l.pos = p;
            if (at_end(&l)) break;
            const char c = src.ptr[l.pos];
            if (c == '#' || c == '\n' || c == '\r') {
                // Not a logical line: indentation does not count
                if (c == '#') {
                    const size_t beg = l.pos;
                    while (!at_end(&l) && !is_newline_at(&l, l.pos)) l.pos += 1;
                    emit(&l, MD_VIAL_TOKEN_COMMENT, beg, l.pos, 0, 0);
                }
                l.pos += newline_len(&l, l.pos);
                continue;
            }
            handle_indent(&l, ind, l.pos);
            at_line_start = false;
        }

        const int c = peek(&l, 0);
        if (c == ' ' || c == '\t' || c == '\f') { l.pos += 1; continue; }

        if (c == '\n' || c == '\r') {
            const size_t beg = l.pos;
            const size_t n = newline_len(&l, l.pos);
            l.pos += n;
            if (md_array_size(l.open) > 0) {
                const uint32_t col = l.open[0].line_col8;
                if (!next_line_is_statement(&l, l.pos, col)) continue;   // a line break inside brackets
                close_all_open(&l);
            }
            emit_newline(&l, beg, beg + n);
            at_line_start = true;
            continue;
        }

        if (c == '#') {
            const size_t beg = l.pos;
            while (!at_end(&l) && !is_newline_at(&l, l.pos)) l.pos += 1;
            emit(&l, MD_VIAL_TOKEN_COMMENT, beg, l.pos, 0, 0);
            continue;
        }

        if (c == '\\') {
            const size_t n = newline_len(&l, l.pos + 1);
            if (n) { l.pos += 1 + n; continue; }   // an explicit line join
            l.pos += 1;
            vial_diag(diags, "E0101", span(l.pos - 1, l.pos), "unexpected `\\`; only a line break may follow it");
            emit(&l, MD_VIAL_TOKEN_ERROR, l.pos - 1, l.pos, 0, 0);
            continue;
        }

        if (is_ident_start(c)) {
            const size_t beg = l.pos;
            while (!at_end(&l) && is_ident_char(peek(&l, 0))) l.pos += 1;
            const str_t word = {src.ptr + beg, l.pos - beg};
            uint16_t flags;
            if (!at_end(&l) && (peek(&l, 0) == '\'' || peek(&l, 0) == '"') && string_prefix(word, &flags)) {
                lex_string(&l, beg, flags);
                continue;
            }
            const vial_kw_t kw = keyword_of(word);
            // A unit name right after a number is the number's unit, whatever the script binds that name to
            const size_t nt = md_array_size(l.tokens);
            const md_vial_token_kind_t prev = nt ? l.tokens[nt - 1].kind : MD_VIAL_TOKEN_EOF;
            if (!kw && (prev == MD_VIAL_TOKEN_INT || prev == MD_VIAL_TOKEN_FLOAT) && vial_unit_name(word, NULL)) {
                emit(&l, MD_VIAL_TOKEN_UNIT, beg, l.pos, 0, 0);
                continue;
            }
            emit(&l, kw ? MD_VIAL_TOKEN_KEYWORD : MD_VIAL_TOKEN_NAME, beg, l.pos, (uint16_t)kw, 0);
            continue;
        }

        if (is_digit(c) || (c == '.' && is_digit(peek(&l, 1)))) {
            lex_number(&l);
            continue;
        }

        if (c == '\'' || c == '"') {
            lex_string(&l, l.pos, 0);
            continue;
        }

        lex_op(&l);
    }

    if (md_array_size(l.open) > 0) close_all_open(&l);
    emit_newline(&l, src.len, src.len);
    while (md_array_size(l.indents) > 1) {
        md_array_pop(l.indents);
        emit(&l, MD_VIAL_TOKEN_DEDENT, src.len, src.len, 0, 0);
    }
    emit(&l, MD_VIAL_TOKEN_EOF, src.len, src.len, 0, 0);

    md_array_free(l.open, alloc);
    md_array_free(l.indents, alloc);
    return l.tokens;
}

// Appends formatted text to an md_array of chars
static void buf_fmt(char** buf, md_allocator_i* alloc, const char* fmt, ...) {
    char tmp[256];
    va_list args;
    va_start(args, fmt);
    int len = vsnprintf(tmp, sizeof(tmp), fmt, args);
    va_end(args);
    if (len < 0) return;
    if ((size_t)len >= sizeof(tmp)) len = (int)sizeof(tmp) - 1;
    md_array_push_array(*buf, tmp, (size_t)len, alloc);
}

str_t vial_tokens_dump(const vial_token_t* tokens, str_t src, md_allocator_i* alloc) {
    char* buf = NULL;
    for (size_t i = 0; i < md_array_size(tokens); ++i) {
        const vial_token_t* t = &tokens[i];
        buf_fmt(&buf, alloc, "%s %u..%u", vial_token_kind_str(t->kind), t->span.beg, t->span.end);
        switch (t->kind) {
        case MD_VIAL_TOKEN_NEWLINE: case MD_VIAL_TOKEN_INDENT: case MD_VIAL_TOKEN_DEDENT: case MD_VIAL_TOKEN_EOF:
            break;
        default:
            buf_fmt(&buf, alloc, " %.*s", (int)(t->span.end - t->span.beg), src.ptr + t->span.beg);
            break;
        }
        buf_fmt(&buf, alloc, "\n");
    }
    return (str_t){buf, md_array_size(buf)};
}

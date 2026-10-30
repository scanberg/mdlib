// The VIAL parser: recursive descent for statements, one function per precedence level for expressions.
//
// It reads the full Python subset of the reference, and more: every construct VIAL refuses is parsed far enough to
// be reported where it is written and skipped cleanly, so no script misreads as something else. Errors never stop
// it. A statement that fails gives one diagnostic, the parser skips to the end of the statement, and the next one
// parses as if nothing happened; an error the lexer already reported (an unclosed bracket, a character that is not
// VIAL) is not reported again. Nesting is bounded (E0104), so no input runs the parser, or a later stage walking the
// tree, out of stack.
//
// Precedence, loosest first, as Python's with `..` added between the comparisons and `|`:
//
//     lambda   if-else   or   and   not   comparisons (< > == >= <= != in, not in, is, is not)
//     ..   |   ^   &   << >>   + -   * / // % @   unary + - ~   **   await   calls, subscripts, attributes

#include "ast.h"

#include <core/md_common.h>

#include <string.h>
#include <stdio.h>
#include <stdarg.h>

typedef struct parser_t {
    vial_token_t*   toks;       // the tokens without comments
    size_t          num;
    size_t          pos;
    uint32_t        last_end;   // where the last consumed token ended
    str_t           src;
    vial_diags_t*   diags;
    md_allocator_i* alloc;
    vial_tree_t     tree;
    vial_node_id*   stack;      // children of the nodes being built
    bool            panic;      // this statement has reported an error; report no more until it ends
    uint32_t        nesting;    // brackets, blocks and prefix operators open around the current token
    uint16_t*       depth;      // md_array, by node id: the height of the node's subtree
} parser_t;

// How deep a script may nest. The parser recurses once per level of brackets, blocks and prefix operators, so
// nesting is bounded by what a small thread stack allows; every later stage walks the tree recursively, so its
// height is bounded too. A long chain such as a + b + c + ... is built in a loop but is as high as it is long.
// Both are far beyond what a script needs, and reported as E0104.
#define VIAL_MAX_NESTING 100
#define VIAL_MAX_HEIGHT  1000

// ### TOKENS ###

static inline const vial_token_t* cur(const parser_t* p) { return &p->toks[p->pos]; }
static inline const vial_token_t* ahead(const parser_t* p, size_t k) {
    const size_t i = p->pos + k;
    return &p->toks[i < p->num ? i : p->num - 1];
}

static inline bool tok_is_op(const vial_token_t* t, vial_op_t op) {
    return (t->kind == MD_VIAL_TOKEN_OPERATOR || t->kind == MD_VIAL_TOKEN_PUNCT) && t->sub == op;
}
static inline bool tok_is_kw(const vial_token_t* t, vial_kw_t kw) { return t->kind == MD_VIAL_TOKEN_KEYWORD && t->sub == kw; }
static inline bool is_op(const parser_t* p, vial_op_t op) { return tok_is_op(cur(p), op); }
static inline bool is_kw(const parser_t* p, vial_kw_t kw) { return tok_is_kw(cur(p), kw); }
static inline bool is_kind(const parser_t* p, md_vial_token_kind_t k) { return cur(p)->kind == k; }
static inline bool at_line_end(const parser_t* p) {
    const md_vial_token_kind_t k = cur(p)->kind;
    return k == MD_VIAL_TOKEN_NEWLINE || k == MD_VIAL_TOKEN_EOF || k == MD_VIAL_TOKEN_DEDENT;
}

static inline void advance(parser_t* p) {
    if (p->pos + 1 < p->num) {
        // Spans end at the last real token: a statement does not take in its line break, a block not its dedent
        const md_vial_token_kind_t k = cur(p)->kind;
        if (k != MD_VIAL_TOKEN_NEWLINE && k != MD_VIAL_TOKEN_INDENT && k != MD_VIAL_TOKEN_DEDENT) p->last_end = cur(p)->span.end;
        p->pos += 1;
    }
}

static inline bool accept_op(parser_t* p, vial_op_t op) { if (is_op(p, op)) { advance(p); return true; } return false; }
static inline bool accept_kw(parser_t* p, vial_kw_t kw) { if (is_kw(p, kw)) { advance(p); return true; } return false; }

static inline md_vial_span_t tok_span(const vial_token_t* t) { return t->span; }

// How a token reads in a message
static void describe(char* buf, size_t cap, const parser_t* p, const vial_token_t* t) {
    switch (t->kind) {
    case MD_VIAL_TOKEN_NEWLINE: snprintf(buf, cap, "end of line"); break;
    case MD_VIAL_TOKEN_EOF:     snprintf(buf, cap, "end of text"); break;
    case MD_VIAL_TOKEN_INDENT:  snprintf(buf, cap, "indent"); break;
    case MD_VIAL_TOKEN_DEDENT:  snprintf(buf, cap, "end of block"); break;
    default: {
        const uint32_t len = t->span.end - t->span.beg;
        snprintf(buf, cap, "`%.*s`", (int)(len > 24 ? 24 : len), p->src.ptr + t->span.beg);
        break;
    }
    }
}

// ### DIAGNOSTICS ###

// Has a diagnostic been given at this position already? The lexer's, for a bracket it closed, or a token it refused.
static bool reported_at(const parser_t* p, uint32_t beg) {
    return vial_diag_at(p->diags, beg);
}

static void error(parser_t* p, const char code[6], md_vial_span_t span, const char* fmt, ...) {
    if (p->panic) return;
    char buf[400];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    vial_diag(p->diags, code, span, "%s", buf);
    p->panic = true;
}

// A refused construct: reported, but the statement goes on parsing
static void refuse(parser_t* p, const char code[6], md_vial_span_t span, const char* fmt, ...) {
    if (p->panic) return;
    char buf[400];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    vial_diag(p->diags, code, span, "%s", buf);
}

// The current token is not what the grammar wants here
static void unexpected(parser_t* p, const char* expected) {
    const vial_token_t* t = cur(p);
    if (t->kind == MD_VIAL_TOKEN_ERROR || (t->flags & VIAL_TOK_REPORTED) || reported_at(p, t->span.beg)) {
        p->panic = true;   // the lexer said so already
        return;
    }
    char what[64];
    describe(what, sizeof(what), p, t);
    if (expected) error(p, "E0101", tok_span(t), "unexpected %s; expected %s", what, expected);
    else          error(p, "E0101", tok_span(t), "unexpected %s", what);
}

static bool expect_op(parser_t* p, vial_op_t op, const char* expected) {
    if (accept_op(p, op)) return true;
    unexpected(p, expected);
    return false;
}

static void expect_close(parser_t* p, const vial_token_t* opener, vial_op_t closer) {
    const vial_token_t* t = cur(p);
    if (tok_is_op(t, closer)) { advance(p); return; }
    // A closer the lexer reported as mismatched still closes
    if (t->kind == MD_VIAL_TOKEN_PUNCT && (t->flags & VIAL_TOK_REPORTED) &&
        (t->sub == VIAL_OP_RPAREN || t->sub == VIAL_OP_RBRACKET || t->sub == VIAL_OP_RBRACE)) {
        advance(p);
        return;
    }
    if (reported_at(p, opener->span.beg)) { p->panic = true; return; }   // never closed: the lexer said so
    char want[48];
    snprintf(want, sizeof(want), "`%s` to close `%s` on line %u", vial_op_str(closer), vial_op_str((vial_op_t)opener->sub),
        vial_line_of(p->diags, opener->span.beg));
    unexpected(p, want);
}

static const char* unsupported_message(vial_unsupported_t u) {
    switch (u) {
    case VIAL_UNSUP_CLASS:       return "classes are not part of VIAL";
    case VIAL_UNSUP_IMPORT:      return "`import` is not part of VIAL; every function is built in or written in the script";
    case VIAL_UNSUP_LAMBDA:      return "`lambda` is not part of VIAL; a nested def does what it would";
    case VIAL_UNSUP_EXCEPTION:   return "exceptions are not part of VIAL";
    case VIAL_UNSUP_GENERATOR:   return "generators are not part of VIAL; a list comprehension builds the array";
    case VIAL_UNSUP_GLOBAL:      return "`global` and `nonlocal` are not part of VIAL";
    case VIAL_UNSUP_DECORATOR:   return "decorators are not part of VIAL";
    case VIAL_UNSUP_WITH:        return "`with` is not part of VIAL";
    case VIAL_UNSUP_ASYNC:       return "`async` and `await` are not part of VIAL";
    case VIAL_UNSUP_DICT:        return "dicts are not part of VIAL";
    case VIAL_UNSUP_SET:         return "set literals are not part of VIAL; sets of particles come from selectors";
    case VIAL_UNSUP_STARRED:     return "`*args` and `**kwargs` are not part of VIAL";
    case VIAL_UNSUP_FSTRING:     return "f-strings are not part of VIAL";
    case VIAL_UNSUP_DEL:         return "`del` is not part of VIAL";
    case VIAL_UNSUP_WALRUS:      return "`:=` is not part of VIAL";
    case VIAL_UNSUP_LOOP_ELSE:   return "`else` after a loop is not part of VIAL";
    case VIAL_UNSUP_ANNOTATION:  return "variable annotations are not part of VIAL; a parameter may take a kind";
    case VIAL_UNSUP_BYTES:       return "bytes literals are not part of VIAL";
    case VIAL_UNSUP_ITEM_ASSIGN: return "item assignment is not supported yet";
    case VIAL_UNSUP_V1_RANGE:    return "a v1 range";
    default:                     return "not part of VIAL";
    }
}

// ### NODES ###

static inline uint32_t mark(const parser_t* p) { return (uint32_t)md_array_size(p->stack); }
static inline void push(parser_t* p, vial_node_id id) { md_array_push(p->stack, id, p->alloc); }

// A node whose children are the ids pushed since m
static vial_node_id make(parser_t* p, vial_node_kind_t kind, uint16_t sub, uint32_t beg, uint32_t end, uint32_t m) {
    vial_node_t n = {.kind = (uint8_t)kind, .sub = sub, .span = {beg, end < beg ? beg : end}};
    const uint32_t count = mark(p) - m;
    // A node covers its children, also where recovery left one beyond the tokens the node consumed (an ERROR node
    // sits where the missing part was wanted, at the next token)
    uint16_t height = 1;
    for (uint32_t i = 0; i < count; ++i) {
        const vial_node_id c = p->stack[m + i];
        if (!c) continue;
        const md_vial_span_t s = p->tree.nodes[c].span;
        if (s.beg < n.span.beg) n.span.beg = s.beg;
        if (s.end > n.span.end) n.span.end = s.end;
        if (p->depth[c] >= height) height = p->depth[c] + 1;
    }
    if (height > VIAL_MAX_HEIGHT) {
        // Reported where the height is crossed, once; the panic stops the loop that is building the chain. Above
        // that the height is held at the limit, so the nodes around it are not reported again.
        bool crossed = true;
        for (uint32_t i = 0; i < count; ++i) if (p->stack[m + i] && p->depth[p->stack[m + i]] > VIAL_MAX_HEIGHT) crossed = false;
        if (crossed && !p->panic) {
            vial_diag(p->diags, "E0104", n.span, "this expression is more than %d operations deep", VIAL_MAX_HEIGHT);
            p->panic = true;
        }
        height = VIAL_MAX_HEIGHT + 1;
    }
    n.child_beg = (uint32_t)md_array_size(p->tree.children);
    n.child_count = count;
    if (count) md_array_push_array(p->tree.children, p->stack + m, count, p->alloc);
    md_array_shrink(p->stack, m);
    const vial_node_id id = (vial_node_id)md_array_size(p->tree.nodes);
    md_array_push(p->tree.nodes, n, p->alloc);
    md_array_push(p->depth, height, p->alloc);
    return id;
}

static vial_node_id leaf(parser_t* p, vial_node_kind_t kind, uint16_t sub, const vial_token_t* t) {
    return make(p, kind, sub, t->span.beg, t->span.end, mark(p));
}

static inline uint32_t node_beg(const parser_t* p, vial_node_id id) { return p->tree.nodes[id].span.beg; }
static inline uint32_t node_end(const parser_t* p, vial_node_id id) { return p->tree.nodes[id].span.end; }
static inline vial_node_kind_t node_kind(const parser_t* p, vial_node_id id) { return (vial_node_kind_t)p->tree.nodes[id].kind; }

static vial_node_id error_node(parser_t* p) {
    const vial_token_t* t = cur(p);
    return make(p, VIAL_N_ERROR, 0, t->span.beg, t->span.beg, mark(p));
}

// Entering one more level of nesting; false, after reporting, when that is one too many
static bool enter(parser_t* p) {
    if (p->nesting >= VIAL_MAX_NESTING) {
        error(p, "E0104", tok_span(cur(p)), "this nests more than %d levels deep", VIAL_MAX_NESTING);
        p->panic = true;
        return false;
    }
    p->nesting += 1;
    return true;
}
static inline void leave(parser_t* p) { p->nesting -= 1; }

// ### EXPRESSIONS ###

static vial_node_id parse_test(parser_t* p);
static vial_node_id parse_test_nested(parser_t* p);
static vial_node_id parse_or(parser_t* p);
static vial_node_id parse_bitor(parser_t* p);
static vial_node_id parse_factor(parser_t* p);
static vial_node_id parse_testlist(parser_t* p);

// Skips to the token that closes the bracket just consumed, for constructs that are refused whole
static void skip_bracketed(parser_t* p) {
    int depth = 1;
    while (!is_kind(p, MD_VIAL_TOKEN_EOF) && !is_kind(p, MD_VIAL_TOKEN_NEWLINE)) {
        const vial_token_t* t = cur(p);
        if (t->kind == MD_VIAL_TOKEN_PUNCT) {
            if (t->sub == VIAL_OP_LPAREN || t->sub == VIAL_OP_LBRACKET || t->sub == VIAL_OP_LBRACE) depth += 1;
            if (t->sub == VIAL_OP_RPAREN || t->sub == VIAL_OP_RBRACKET || t->sub == VIAL_OP_RBRACE) {
                depth -= 1;
                if (depth == 0) { advance(p); return; }
            }
        }
        advance(p);
    }
}

// Comprehension clauses after the element: for target in iterable [if condition]... [for ...]
static void parse_comp_fors(parser_t* p) {
    while (is_kw(p, VIAL_KW_FOR) || is_kw(p, VIAL_KW_ASYNC)) {
        const uint32_t beg = cur(p)->span.beg;
        if (is_kw(p, VIAL_KW_ASYNC)) { refuse(p, "E0803", tok_span(cur(p)), "%s", unsupported_message(VIAL_UNSUP_ASYNC)); advance(p); }
        if (!accept_kw(p, VIAL_KW_FOR)) { unexpected(p, "`for`"); return; }
        const uint32_t m = mark(p);
        // The target: names and tuples of them, read at the level of | so that `in` ends it
        {
            const uint32_t tm = mark(p);
            const uint32_t tbeg = cur(p)->span.beg;
            push(p, parse_bitor(p));
            bool tuple = false;
            while (accept_op(p, VIAL_OP_COMMA)) {
                tuple = true;
                if (is_kw(p, VIAL_KW_IN)) break;
                push(p, parse_bitor(p));
            }
            if (tuple) push(p, make(p, VIAL_N_TUPLE, 0, tbeg, p->last_end, tm));
        }
        if (!accept_kw(p, VIAL_KW_IN)) { unexpected(p, "`in`"); push(p, error_node(p)); }
        else push(p, parse_or(p));
        while (accept_kw(p, VIAL_KW_IF)) push(p, parse_or(p));
        push(p, make(p, VIAL_N_COMP_FOR, 0, beg, p->last_end, m));
    }
}

// Adjacent string literals are one string
static vial_node_id parse_strings(parser_t* p) {
    const vial_token_t* first = cur(p);
    uint16_t flags = 0;
    while (is_kind(p, MD_VIAL_TOKEN_STRING)) {
        const vial_token_t* t = cur(p);
        if ((t->flags & VIAL_STR_FORMAT) && !(flags & VIAL_STR_FORMAT)) refuse(p, "E0803", tok_span(t), "%s", unsupported_message(VIAL_UNSUP_FSTRING));
        if ((t->flags & VIAL_STR_BYTES) && !(flags & VIAL_STR_BYTES)) refuse(p, "E0803", tok_span(t), "%s", unsupported_message(VIAL_UNSUP_BYTES));
        flags |= t->flags;
        advance(p);
    }
    vial_node_id id = make(p, VIAL_N_STRING, 0, first->span.beg, p->last_end, mark(p));
    p->tree.nodes[id].flags = (uint8_t)(flags & 0xFF);
    return id;
}

// ( ... ): a parenthesized expression, a tuple, or a generator
static vial_node_id parse_paren(parser_t* p) {
    const vial_token_t* open = cur(p);
    const uint32_t beg = open->span.beg;
    advance(p);
    const uint32_t m = mark(p);
    if (is_op(p, VIAL_OP_RPAREN)) {
        advance(p);
        return make(p, VIAL_N_TUPLE, 0, beg, p->last_end, m);
    }
    vial_node_id first = parse_test(p);
    if (is_kw(p, VIAL_KW_FOR)) {
        refuse(p, "E0803", (md_vial_span_t){beg, node_end(p, first)}, "%s", unsupported_message(VIAL_UNSUP_GENERATOR));
        push(p, first);
        parse_comp_fors(p);
        expect_close(p, open, VIAL_OP_RPAREN);
        vial_node_id comp = make(p, VIAL_N_LISTCOMP, 0, beg, p->last_end, m);
        const uint32_t m3 = mark(p);
        push(p, comp);
        return make(p, VIAL_N_UNSUPPORTED, VIAL_UNSUP_GENERATOR, beg, p->last_end, m3);
    }
    if (!is_op(p, VIAL_OP_COMMA)) {
        expect_close(p, open, VIAL_OP_RPAREN);
        return first;   // parentheses only group
    }
    push(p, first);
    while (accept_op(p, VIAL_OP_COMMA)) {
        if (is_op(p, VIAL_OP_RPAREN)) break;
        push(p, parse_test(p));
        if (p->panic) break;
    }
    expect_close(p, open, VIAL_OP_RPAREN);
    return make(p, VIAL_N_TUPLE, 0, beg, p->last_end, m);
}

// [ ... ]: a list, or a list comprehension
static vial_node_id parse_list(parser_t* p) {
    const vial_token_t* open = cur(p);
    const uint32_t beg = open->span.beg;
    advance(p);
    const uint32_t m = mark(p);
    if (is_op(p, VIAL_OP_RBRACKET)) {
        advance(p);
        return make(p, VIAL_N_LIST, 0, beg, p->last_end, m);
    }
    push(p, parse_test(p));
    if (is_kw(p, VIAL_KW_FOR) || is_kw(p, VIAL_KW_ASYNC)) {
        parse_comp_fors(p);
        expect_close(p, open, VIAL_OP_RBRACKET);
        return make(p, VIAL_N_LISTCOMP, 0, beg, p->last_end, m);
    }
    while (accept_op(p, VIAL_OP_COMMA)) {
        if (is_op(p, VIAL_OP_RBRACKET)) break;
        push(p, parse_test(p));
        if (p->panic) break;
    }
    expect_close(p, open, VIAL_OP_RBRACKET);
    return make(p, VIAL_N_LIST, 0, beg, p->last_end, m);
}

// { ... }: a dict or a set, both refused, and skipped whole
static vial_node_id parse_brace(parser_t* p) {
    const uint32_t beg = cur(p)->span.beg;
    advance(p);
    // A dict when it is empty or a ':' follows its first element at depth 1
    bool dict = is_op(p, VIAL_OP_RBRACE);
    for (size_t i = p->pos, depth = 1; i < p->num && !dict; ++i) {
        const vial_token_t* t = &p->toks[i];
        if (t->kind == MD_VIAL_TOKEN_NEWLINE || t->kind == MD_VIAL_TOKEN_EOF) break;
        if (t->kind == MD_VIAL_TOKEN_PUNCT) {
            if (t->sub == VIAL_OP_LPAREN || t->sub == VIAL_OP_LBRACKET || t->sub == VIAL_OP_LBRACE) depth += 1;
            else if (t->sub == VIAL_OP_RPAREN || t->sub == VIAL_OP_RBRACKET || t->sub == VIAL_OP_RBRACE) { if (--depth == 0) break; }
            else if (depth == 1 && t->sub == VIAL_OP_COLON) dict = true;
            else if (depth == 1 && t->sub == VIAL_OP_COMMA) break;
        }
    }
    skip_bracketed(p);
    const vial_unsupported_t u = dict ? VIAL_UNSUP_DICT : VIAL_UNSUP_SET;
    refuse(p, "E0803", (md_vial_span_t){beg, p->last_end}, "%s", unsupported_message(u));
    return make(p, VIAL_N_UNSUPPORTED, (uint16_t)u, beg, p->last_end, mark(p));
}

static vial_node_id parse_atom_(parser_t* p);

// Every bracket, and every prefix that VIAL refuses, recurses through here: one level of nesting each
static vial_node_id parse_atom(parser_t* p) {
    if (!enter(p)) return error_node(p);
    const vial_node_id e = parse_atom_(p);
    leave(p);
    return e;
}

// A number, and its unit if one follows: a name from the unit table, `10 ns`, or a string, `1.5 'kJ/mol'`. The
// quantity is one atom, so t / 2 ns is t / (2 ns).
static vial_node_id parse_number(parser_t* p) {
    const vial_token_t* t = cur(p);
    advance(p);
    const vial_node_id num = leaf(p, t->kind == MD_VIAL_TOKEN_INT ? VIAL_N_INT : VIAL_N_FLOAT, 0, t);
    const vial_token_t* u = cur(p);
    if (u->kind == MD_VIAL_TOKEN_UNIT) {
        advance(p);
        const uint32_t m = mark(p);
        push(p, num);
        push(p, leaf(p, VIAL_N_UNIT, 0, u));
        return make(p, VIAL_N_QUANTITY, 0, t->span.beg, u->span.end, m);
    }
    if (u->kind == MD_VIAL_TOKEN_STRING && !(u->flags & (VIAL_STR_FORMAT | VIAL_STR_BYTES))) {
        advance(p);
        if (!(u->flags & VIAL_STR_UNCLOSED)) {
            // The text between the quotes
            const uint32_t q = (u->flags & VIAL_STR_TRIPLE) ? 3 : 1;
            uint32_t b = u->span.beg;
            while (b < u->span.end && p->src.ptr[b] != '\'' && p->src.ptr[b] != '"') ++b;   // past r
            const str_t text = {p->src.ptr + b + q, (u->span.end - b) >= 2 * q ? u->span.end - b - 2 * q : 0};
            if (!vial_unit_string(text, NULL)) {
                refuse(p, "E0502", tok_span(u), "`%.*s` is not a unit", (int)(u->span.end - u->span.beg), p->src.ptr + u->span.beg);
            }
        }
        const uint32_t m = mark(p);
        push(p, num);
        push(p, leaf(p, VIAL_N_UNIT, 1, u));
        return make(p, VIAL_N_QUANTITY, 0, t->span.beg, u->span.end, m);
    }
    return num;
}

static vial_node_id parse_atom_(parser_t* p) {
    const vial_token_t* t = cur(p);
    switch (t->kind) {
    case MD_VIAL_TOKEN_NAME:   advance(p); return leaf(p, VIAL_N_NAME, 0, t);
    case MD_VIAL_TOKEN_INT:
    case MD_VIAL_TOKEN_FLOAT:  return parse_number(p);
    case MD_VIAL_TOKEN_STRING: return parse_strings(p);
    case MD_VIAL_TOKEN_ERROR:  advance(p); p->panic = true; return leaf(p, VIAL_N_ERROR, 0, t);  // reported by the lexer
    case MD_VIAL_TOKEN_KEYWORD:
        if (t->sub == VIAL_KW_TRUE || t->sub == VIAL_KW_FALSE || t->sub == VIAL_KW_NONE_) {
            advance(p);
            return leaf(p, VIAL_N_CONST, t->sub, t);
        }
        if (t->sub == VIAL_KW_YIELD || t->sub == VIAL_KW_AWAIT) {
            const vial_unsupported_t u = t->sub == VIAL_KW_YIELD ? VIAL_UNSUP_GENERATOR : VIAL_UNSUP_ASYNC;
            refuse(p, "E0803", tok_span(t), "%s", unsupported_message(u));
            advance(p);
            const uint32_t m = mark(p);
            if (!at_line_end(p) && !is_op(p, VIAL_OP_RPAREN)) push(p, u == VIAL_UNSUP_ASYNC ? parse_factor(p) : parse_testlist(p));
            return make(p, VIAL_N_UNSUPPORTED, (uint16_t)u, t->span.beg, p->last_end, m);
        }
        break;
    case MD_VIAL_TOKEN_OPERATOR:
    case MD_VIAL_TOKEN_PUNCT:
        switch (t->sub) {
        case VIAL_OP_LPAREN:   return parse_paren(p);
        case VIAL_OP_LBRACKET: return parse_list(p);
        case VIAL_OP_LBRACE:   return parse_brace(p);
        case VIAL_OP_ELLIPSIS: advance(p); return leaf(p, VIAL_N_ELLIPSIS, 0, t);
        case VIAL_OP_STAR:
        case VIAL_OP_POW: {
            // *x is refused; a '*' with nothing after it is just misplaced
            const vial_token_t* n = ahead(p, 1);
            if (n->kind == MD_VIAL_TOKEN_NEWLINE || n->kind == MD_VIAL_TOKEN_EOF || n->kind == MD_VIAL_TOKEN_PUNCT) break;
            refuse(p, "E0803", tok_span(t), "%s", unsupported_message(VIAL_UNSUP_STARRED));
            advance(p);
            const uint32_t m = mark(p);
            push(p, parse_bitor(p));
            return make(p, VIAL_N_UNSUPPORTED, VIAL_UNSUP_STARRED, t->span.beg, p->last_end, m);
        }
        default: break;
        }
        break;
    default: break;
    }
    unexpected(p, "an expression");
    return error_node(p);
}

// The arguments of a call, after its '('. Returns whether a v1 range such as residue(1:4) was among them.
static bool parse_call_args(parser_t* p, const vial_token_t* open) {
    bool v1 = false;
    while (!is_op(p, VIAL_OP_RPAREN) && !at_line_end(p) && !p->panic) {
        const vial_token_t* t = cur(p);
        if (tok_is_op(t, VIAL_OP_STAR) || tok_is_op(t, VIAL_OP_POW)) {
            push(p, parse_atom(p));   // refused there
        } else if (t->kind == MD_VIAL_TOKEN_NAME && tok_is_op(ahead(p, 1), VIAL_OP_ASSIGN)) {
            const uint32_t m = mark(p);
            advance(p);
            push(p, leaf(p, VIAL_N_NAME, 0, t));
            advance(p);   // =
            push(p, parse_test(p));
            push(p, make(p, VIAL_N_KWARG, 0, t->span.beg, p->last_end, m));
        } else if (tok_is_op(t, VIAL_OP_COLON)) {
            // v1: residue(:), residue(:4)
            const uint32_t m = mark(p);
            push(p, 0);
            advance(p);
            push(p, (is_op(p, VIAL_OP_COMMA) || is_op(p, VIAL_OP_RPAREN)) ? 0 : parse_test(p));
            push(p, make(p, VIAL_N_UNSUPPORTED, VIAL_UNSUP_V1_RANGE, t->span.beg, p->last_end, m));
            v1 = true;
        } else {
            vial_node_id arg = parse_test(p);
            if (is_kw(p, VIAL_KW_FOR)) {
                refuse(p, "E0803", (md_vial_span_t){node_beg(p, arg), node_end(p, arg)}, "%s", unsupported_message(VIAL_UNSUP_GENERATOR));
                const uint32_t m = mark(p);
                push(p, arg);
                parse_comp_fors(p);
                vial_node_id comp = make(p, VIAL_N_LISTCOMP, 0, node_beg(p, arg), p->last_end, m);
                const uint32_t m2 = mark(p);
                push(p, comp);
                arg = make(p, VIAL_N_UNSUPPORTED, VIAL_UNSUP_GENERATOR, node_beg(p, arg), p->last_end, m2);
            } else if (is_op(p, VIAL_OP_COLON)) {
                // v1: residue(1:4)
                const uint32_t m = mark(p);
                push(p, arg);
                advance(p);
                push(p, (is_op(p, VIAL_OP_COMMA) || is_op(p, VIAL_OP_RPAREN)) ? 0 : parse_test(p));
                arg = make(p, VIAL_N_UNSUPPORTED, VIAL_UNSUP_V1_RANGE, node_beg(p, arg), p->last_end, m);
                v1 = true;
            }
            push(p, arg);
        }
        if (!accept_op(p, VIAL_OP_COMMA)) break;
    }
    expect_close(p, open, VIAL_OP_RPAREN);
    return v1;
}

// One subscript: an expression or a slice a:b:c with any part left out
static vial_node_id parse_subscript(parser_t* p) {
    const uint32_t beg = cur(p)->span.beg;
    vial_node_id start = 0;
    if (!is_op(p, VIAL_OP_COLON)) {
        start = parse_test(p);
        if (!is_op(p, VIAL_OP_COLON)) return start;
    }
    const uint32_t m = mark(p);
    push(p, start);
    advance(p);   // :
    const bool stop_absent = is_op(p, VIAL_OP_COLON) || is_op(p, VIAL_OP_COMMA) || is_op(p, VIAL_OP_RBRACKET);
    push(p, stop_absent ? 0 : parse_test(p));
    if (accept_op(p, VIAL_OP_COLON)) {
        const bool step_absent = is_op(p, VIAL_OP_COMMA) || is_op(p, VIAL_OP_RBRACKET);
        push(p, step_absent ? 0 : parse_test(p));
    } else {
        push(p, 0);
    }
    return make(p, VIAL_N_SLICE, 0, beg, p->last_end, m);
}

// v1's selectors by 1-based inclusive position, and what VIAL writes for them
static const char* v1_plural(str_t name) {
    static const char* map[][2] = {
        {"residue", "residues"}, {"chain", "chains"}, {"atom", "particles"}, {"molecule", "molecules"},
        {"part", "parts"}, {"polymer", "parts"}, {"protein", "parts"}, {"nucleic", "parts"},
    };
    for (size_t i = 0; i < ARRAY_SIZE(map); ++i) if (str_eq_cstr(name, map[i][0])) return map[i][1];
    return NULL;
}

static bool int_literal(const parser_t* p, vial_node_id id, int64_t* out) {
    if (!id || node_kind(p, id) != VIAL_N_INT) return false;
    str_t s = {p->src.ptr + node_beg(p, id), node_end(p, id) - node_beg(p, id)};
    int64_t v = 0;
    for (size_t i = 0; i < s.len; ++i) {
        if (s.ptr[i] == '_') continue;
        if (!is_digit(s.ptr[i]) || v > 1000000000) return false;
        v = v * 10 + (s.ptr[i] - '0');
    }
    *out = v;
    return true;
}

// E0701 on a v1 call such as residue(1:4), with the VIAL spelling as the fix where it is clear: residues[0:4]
static void report_v1_call(parser_t* p, vial_node_id call) {
    const md_vial_span_t s = p->tree.nodes[call].span;
    const uint32_t len = s.end - s.beg;
    const size_t before = md_array_size(p->diags->items);
    refuse(p, "E0701", s, "`%.*s%s` is v1; positions are subscripts", (int)(len > 40 ? 40 : len), p->src.ptr + s.beg,
        len > 40 ? "…" : "");
    if (md_array_size(p->diags->items) == before) return;

    // The fix: one argument, a v1 range of int literals or open ends, on a selector v1 had
    const vial_node_t* n = &p->tree.nodes[call];
    const vial_node_id callee = vial_child(&p->tree, call, 0);
    if (n->child_count != 2 || node_kind(p, callee) != VIAL_N_NAME) return;
    const vial_node_id arg = vial_child(&p->tree, call, 1);
    if (node_kind(p, arg) != VIAL_N_UNSUPPORTED || p->tree.nodes[arg].sub != VIAL_UNSUP_V1_RANGE) return;
    const char* plural = v1_plural((str_t){p->src.ptr + node_beg(p, callee), node_end(p, callee) - node_beg(p, callee)});
    if (!plural) return;
    const vial_node_id lo = vial_child(&p->tree, arg, 0), hi = vial_child(&p->tree, arg, 1);
    int64_t a = 0, b = 0;
    if ((lo && !int_literal(p, lo, &a)) || (hi && !int_literal(p, hi, &b)) || (lo && a < 1)) return;
    char fix[96];
    if (!lo && !hi)    snprintf(fix, sizeof(fix), "%s", plural);
    else if (!hi)      snprintf(fix, sizeof(fix), "%s[%lld:]", plural, (long long)(a - 1));
    else if (!lo)      snprintf(fix, sizeof(fix), "%s[:%lld]", plural, (long long)b);
    else               snprintf(fix, sizeof(fix), "%s[%lld:%lld]", plural, (long long)(a - 1), (long long)b);
    vial_diag_set_fix(p->diags, before, str_from_cstr(fix));
}

static vial_node_id parse_primary(parser_t* p) {
    vial_node_id e = parse_atom(p);
    for (;;) {
        if (p->panic) return e;   // an atom that failed takes no call, subscript or attribute
        const vial_token_t* t = cur(p);
        if (tok_is_op(t, VIAL_OP_LPAREN)) {
            const uint32_t m = mark(p);
            push(p, e);
            advance(p);
            const bool v1 = parse_call_args(p, t);
            e = make(p, VIAL_N_CALL, 0, node_beg(p, e), p->last_end, m);
            if (v1) report_v1_call(p, e);
        } else if (tok_is_op(t, VIAL_OP_LBRACKET)) {
            const uint32_t m = mark(p);
            push(p, e);
            advance(p);
            const uint32_t ibeg = cur(p)->span.beg;
            const uint32_t im = mark(p);
            push(p, parse_subscript(p));
            bool tuple = false;
            while (accept_op(p, VIAL_OP_COMMA)) {
                tuple = true;
                if (is_op(p, VIAL_OP_RBRACKET)) break;
                push(p, parse_subscript(p));
                if (p->panic) break;
            }
            if (tuple) push(p, make(p, VIAL_N_TUPLE, 0, ibeg, p->last_end, im));
            expect_close(p, t, VIAL_OP_RBRACKET);
            e = make(p, VIAL_N_SUBSCRIPT, 0, node_beg(p, e), p->last_end, m);
        } else if (tok_is_op(t, VIAL_OP_DOT)) {
            advance(p);
            const vial_token_t* name = cur(p);
            const uint32_t m = mark(p);
            push(p, e);
            if (name->kind == MD_VIAL_TOKEN_NAME) {
                advance(p);
                push(p, leaf(p, VIAL_N_NAME, 0, name));
            } else {
                unexpected(p, "a name after `.`");
                push(p, error_node(p));
            }
            e = make(p, VIAL_N_ATTR, 0, node_beg(p, e), p->last_end, m);
        } else {
            return e;
        }
        if (p->panic) return e;
    }
}

static vial_node_id parse_power(parser_t* p) {
    vial_node_id base = parse_primary(p);
    if (is_op(p, VIAL_OP_POW)) {
        advance(p);
        const uint32_t m = mark(p);
        push(p, base);
        if (!enter(p)) {
            push(p, error_node(p));
        } else {
            push(p, parse_factor(p));   // right-associative, and binds tighter than a unary minus on its left
            leave(p);
        }
        return make(p, VIAL_N_BINARY, VIAL_OP_POW, node_beg(p, base), p->last_end, m);
    }
    return base;
}

static vial_node_id parse_factor(parser_t* p) {
    const vial_token_t* t = cur(p);
    if (tok_is_op(t, VIAL_OP_PLUS) || tok_is_op(t, VIAL_OP_MINUS) || tok_is_op(t, VIAL_OP_TILDE)) {
        advance(p);
        if (!enter(p)) return error_node(p);
        const uint32_t m = mark(p);
        push(p, parse_factor(p));
        leave(p);
        return make(p, VIAL_N_UNARY, t->sub, t->span.beg, p->last_end, m);
    }
    return parse_power(p);
}

typedef vial_node_id (*parse_fn)(parser_t*);

// A left-associative level of binary operators
static vial_node_id parse_binary(parser_t* p, parse_fn next, const vial_op_t* ops, size_t num_ops) {
    vial_node_id left = next(p);
    for (;;) {
        const vial_token_t* t = cur(p);
        bool match = false;
        for (size_t i = 0; i < num_ops; ++i) if (tok_is_op(t, ops[i])) { match = true; break; }
        if (!match || p->panic) return left;
        advance(p);
        const uint32_t m = mark(p);
        push(p, left);
        push(p, next(p));
        left = make(p, VIAL_N_BINARY, t->sub, node_beg(p, left), p->last_end, m);
    }
}

static vial_node_id parse_term(parser_t* p) {
    static const vial_op_t ops[] = {VIAL_OP_STAR, VIAL_OP_SLASH, VIAL_OP_DSLASH, VIAL_OP_PERCENT, VIAL_OP_AT};
    return parse_binary(p, parse_factor, ops, ARRAY_SIZE(ops));
}
static vial_node_id parse_arith(parser_t* p) {
    static const vial_op_t ops[] = {VIAL_OP_PLUS, VIAL_OP_MINUS};
    return parse_binary(p, parse_term, ops, ARRAY_SIZE(ops));
}
static vial_node_id parse_shift(parser_t* p) {
    static const vial_op_t ops[] = {VIAL_OP_LSHIFT, VIAL_OP_RSHIFT};
    return parse_binary(p, parse_arith, ops, ARRAY_SIZE(ops));
}
static vial_node_id parse_bitand(parser_t* p) {
    static const vial_op_t ops[] = {VIAL_OP_AMP};
    return parse_binary(p, parse_shift, ops, ARRAY_SIZE(ops));
}
static vial_node_id parse_bitxor(parser_t* p) {
    static const vial_op_t ops[] = {VIAL_OP_CARET};
    return parse_binary(p, parse_bitand, ops, ARRAY_SIZE(ops));
}
static vial_node_id parse_bitor(parser_t* p) {
    static const vial_op_t ops[] = {VIAL_OP_PIPE};
    return parse_binary(p, parse_bitxor, ops, ARRAY_SIZE(ops));
}

// a..b, the inclusive range: looser than arithmetic and bits, tighter than comparisons
static vial_node_id parse_range(parser_t* p) {
    vial_node_id lo = parse_bitor(p);
    if (is_op(p, VIAL_OP_RANGE)) {
        advance(p);
        const uint32_t m = mark(p);
        push(p, lo);
        push(p, parse_bitor(p));
        return make(p, VIAL_N_RANGE, 0, node_beg(p, lo), p->last_end, m);
    }
    return lo;
}

static vial_cmp_t comparison_at(const parser_t* p, size_t* len) {
    const vial_token_t* t = cur(p);
    *len = 1;
    if (t->kind == MD_VIAL_TOKEN_OPERATOR) {
        switch (t->sub) {
        case VIAL_OP_LT: return VIAL_CMP_LT;
        case VIAL_OP_GT: return VIAL_CMP_GT;
        case VIAL_OP_LE: return VIAL_CMP_LE;
        case VIAL_OP_GE: return VIAL_CMP_GE;
        case VIAL_OP_EQ: return VIAL_CMP_EQ;
        case VIAL_OP_NE: return VIAL_CMP_NE;
        default: return 0;
        }
    }
    if (tok_is_kw(t, VIAL_KW_IN)) return VIAL_CMP_IN;
    if (tok_is_kw(t, VIAL_KW_NOT) && tok_is_kw(ahead(p, 1), VIAL_KW_IN)) { *len = 2; return VIAL_CMP_NOT_IN; }
    if (tok_is_kw(t, VIAL_KW_IS)) {
        if (tok_is_kw(ahead(p, 1), VIAL_KW_NOT)) { *len = 2; return VIAL_CMP_IS_NOT; }
        return VIAL_CMP_IS;
    }
    return 0;
}

static vial_node_id parse_comparison(parser_t* p) {
    vial_node_id first = parse_range(p);
    size_t len;
    vial_cmp_t cmp = comparison_at(p, &len);
    if (!cmp) return first;
    const uint32_t m = mark(p);
    push(p, first);
    while (cmp && !p->panic) {
        const vial_token_t* t = cur(p);
        const uint32_t beg = t->span.beg;
        for (size_t i = 0; i < len; ++i) advance(p);
        push(p, make(p, VIAL_N_CMPOP, (uint16_t)cmp, beg, p->last_end, mark(p)));
        push(p, parse_range(p));
        cmp = comparison_at(p, &len);
    }
    return make(p, VIAL_N_COMPARE, 0, node_beg(p, first), p->last_end, m);
}

static vial_node_id parse_not(parser_t* p) {
    const vial_token_t* t = cur(p);
    if (tok_is_kw(t, VIAL_KW_NOT)) {
        advance(p);
        if (!enter(p)) return error_node(p);
        const uint32_t m = mark(p);
        push(p, parse_not(p));
        leave(p);
        return make(p, VIAL_N_NOT, 0, t->span.beg, p->last_end, m);
    }
    return parse_comparison(p);
}

static vial_node_id parse_boolop(parser_t* p, vial_kw_t kw, parse_fn next) {
    vial_node_id first = next(p);
    if (!is_kw(p, kw)) return first;
    const uint32_t m = mark(p);
    push(p, first);
    while (accept_kw(p, kw) && !p->panic) push(p, next(p));
    return make(p, VIAL_N_BOOLOP, (uint16_t)kw, node_beg(p, first), p->last_end, m);
}

static vial_node_id parse_and(parser_t* p) { return parse_boolop(p, VIAL_KW_AND, parse_not); }
static vial_node_id parse_or(parser_t* p)  { return parse_boolop(p, VIAL_KW_OR, parse_and); }

static vial_node_id parse_lambda(parser_t* p) {
    const vial_token_t* t = cur(p);
    refuse(p, "E0803", tok_span(t), "%s", unsupported_message(VIAL_UNSUP_LAMBDA));
    advance(p);
    while (!is_op(p, VIAL_OP_COLON) && !at_line_end(p)) advance(p);   // the parameters
    const uint32_t m = mark(p);
    if (accept_op(p, VIAL_OP_COLON)) push(p, parse_test_nested(p));
    return make(p, VIAL_N_UNSUPPORTED, VIAL_UNSUP_LAMBDA, t->span.beg, p->last_end, m);
}

// parse_test, one level of nesting deeper: for the parts that recurse on the right without passing an atom
static vial_node_id parse_test_nested(parser_t* p) {
    if (!enter(p)) return error_node(p);
    const vial_node_id e = parse_test(p);
    leave(p);
    return e;
}

static vial_node_id parse_test(parser_t* p) {
    if (is_kw(p, VIAL_KW_LAMBDA)) return parse_lambda(p);
    vial_node_id e = parse_or(p);
    if (is_kw(p, VIAL_KW_IF) && !p->panic) {
        advance(p);
        const uint32_t m = mark(p);
        push(p, e);
        push(p, parse_or(p));
        if (accept_kw(p, VIAL_KW_ELSE)) push(p, parse_test_nested(p));
        else { unexpected(p, "`else`"); push(p, error_node(p)); }
        return make(p, VIAL_N_TERNARY, 0, node_beg(p, e), p->last_end, m);
    }
    if (is_op(p, VIAL_OP_WALRUS) && !p->panic) {
        refuse(p, "E0803", tok_span(cur(p)), "%s", unsupported_message(VIAL_UNSUP_WALRUS));
        advance(p);
        const uint32_t m = mark(p);
        push(p, e);
        push(p, parse_test_nested(p));
        return make(p, VIAL_N_UNSUPPORTED, VIAL_UNSUP_WALRUS, node_beg(p, e), p->last_end, m);
    }
    return e;
}

// test, test, ...: a tuple when there is a comma
static vial_node_id parse_testlist(parser_t* p) {
    vial_node_id first = parse_test(p);
    if (!is_op(p, VIAL_OP_COMMA) || p->panic) return first;
    const uint32_t m = mark(p);
    push(p, first);
    while (accept_op(p, VIAL_OP_COMMA)) {
        if (at_line_end(p) || is_op(p, VIAL_OP_ASSIGN) || is_op(p, VIAL_OP_SEMICOLON) || is_op(p, VIAL_OP_RPAREN)) break;
        push(p, parse_test(p));
        if (p->panic) break;
    }
    return make(p, VIAL_N_TUPLE, 0, node_beg(p, first), p->last_end, m);
}

// ### STATEMENTS ###

static void parse_stmt(parser_t* p);
static vial_node_id parse_suite(parser_t* p, const vial_token_t* header);

// After an error: skip to the end of the statement. A block that follows it is parsed, so the statements in it
// are checked, and an error in a header does not turn its whole block into "unexpected indent".
static void synchronize(parser_t* p) {
    while (!is_kind(p, MD_VIAL_TOKEN_NEWLINE) && !is_kind(p, MD_VIAL_TOKEN_EOF) && !is_kind(p, MD_VIAL_TOKEN_DEDENT)) advance(p);
    if (is_kind(p, MD_VIAL_TOKEN_NEWLINE)) advance(p);
    if (is_kind(p, MD_VIAL_TOKEN_INDENT)) {
        p->panic = false;
        push(p, parse_suite(p, NULL));
    }
}

// Statements until the end of a block (a DEDENT, consumed) or of the text
static void parse_stmts_(parser_t* p, bool until_dedent);

// A block nested one too many: skip to its end, nested blocks included
static void skip_block(parser_t* p) {
    int level = 1;
    while (!is_kind(p, MD_VIAL_TOKEN_EOF)) {
        if (is_kind(p, MD_VIAL_TOKEN_INDENT)) level += 1;
        if (is_kind(p, MD_VIAL_TOKEN_DEDENT) && --level == 0) { advance(p); return; }
        advance(p);
    }
}

static void parse_stmts(parser_t* p, bool until_dedent) {
    if (!until_dedent) { parse_stmts_(p, false); return; }
    p->panic = false;
    if (!enter(p)) { skip_block(p); p->panic = false; return; }
    parse_stmts_(p, true);
    leave(p);
}

static void parse_stmts_(parser_t* p, bool until_dedent) {
    for (;;) {
        if (is_kind(p, MD_VIAL_TOKEN_EOF)) return;
        if (is_kind(p, MD_VIAL_TOKEN_DEDENT)) {
            if (until_dedent) { advance(p); return; }
            advance(p);   // cannot happen at the top level, but never loop on it
            continue;
        }
        if (is_kind(p, MD_VIAL_TOKEN_NEWLINE)) { advance(p); continue; }
        if (is_kind(p, MD_VIAL_TOKEN_INDENT)) {
            // An indented line where no block starts: report it, and read its statements as if they were not
            p->panic = false;
            error(p, "E0103", tok_span(cur(p)), "unexpected indent; no block starts on the line before");
            advance(p);
            p->panic = false;
            parse_stmts(p, true);
            continue;
        }
        const size_t before = p->pos;
        parse_stmt(p);
        if (p->pos == before) advance(p);   // never stall
    }
}

static vial_node_id parse_block(parser_t* p, uint32_t beg) {
    const uint32_t m = mark(p);
    parse_stmts(p, true);
    return make(p, VIAL_N_BLOCK, 0, beg, p->last_end, m);
}

// What follows a header's ':' : statements on the same line, or an indented block
static vial_node_id parse_suite(parser_t* p, const vial_token_t* header) {
    const uint32_t beg = cur(p)->span.beg;
    if (is_kind(p, MD_VIAL_TOKEN_NEWLINE)) {
        advance(p);
        if (!is_kind(p, MD_VIAL_TOKEN_INDENT)) {
            char what[64] = "the line before";
            if (header) describe(what, sizeof(what), p, header);
            error(p, "E0103", tok_span(cur(p)), "expected an indented block after %s", what);
            return make(p, VIAL_N_BLOCK, 0, beg, beg, mark(p));
        }
        advance(p);
        p->panic = false;
        return parse_block(p, cur(p)->span.beg);
    }
    if (is_kind(p, MD_VIAL_TOKEN_INDENT)) {   // after synchronize
        advance(p);
        p->panic = false;
        return parse_block(p, cur(p)->span.beg);
    }
    // One line: if a: pass
    const uint32_t m = mark(p);
    if (enter(p)) {
        parse_stmt(p);
        leave(p);
    } else {
        push(p, error_node(p));
        while (!at_line_end(p)) advance(p);   // the rest of the line, or each level above would try it again
    }
    return make(p, VIAL_N_BLOCK, 0, beg, p->last_end, m);
}

// header ':' suite, where the header has been read up to the ':'
static vial_node_id parse_colon_suite(parser_t* p, const vial_token_t* header) {
    if (!expect_op(p, VIAL_OP_COLON, "`:`")) {
        while (!is_op(p, VIAL_OP_COLON) && !at_line_end(p)) advance(p);
        if (!accept_op(p, VIAL_OP_COLON)) {
            if (is_kind(p, MD_VIAL_TOKEN_NEWLINE)) advance(p);
            const uint32_t beg = cur(p)->span.beg;
            if (is_kind(p, MD_VIAL_TOKEN_INDENT)) { advance(p); p->panic = false; return parse_block(p, cur(p)->span.beg); }
            return make(p, VIAL_N_BLOCK, 0, beg, beg, mark(p));
        }
    }
    return parse_suite(p, header);
}

// A refused compound statement: report it, skip its header, parse its body
static void parse_refused_compound(parser_t* p, vial_unsupported_t u) {
    const vial_token_t* t = cur(p);
    refuse(p, "E0803", tok_span(t), "%s", unsupported_message(u));
    advance(p);
    const uint32_t m = mark(p);
    while (!is_op(p, VIAL_OP_COLON) && !at_line_end(p)) advance(p);
    if (accept_op(p, VIAL_OP_COLON)) push(p, parse_suite(p, t));
    // try: also except, else and finally
    while (u == VIAL_UNSUP_EXCEPTION && (is_kw(p, VIAL_KW_EXCEPT) || is_kw(p, VIAL_KW_FINALLY) || is_kw(p, VIAL_KW_ELSE))) {
        const vial_token_t* h = cur(p);
        while (!is_op(p, VIAL_OP_COLON) && !at_line_end(p)) advance(p);
        if (accept_op(p, VIAL_OP_COLON)) push(p, parse_suite(p, h)); else break;
    }
    push(p, make(p, VIAL_N_UNSUPPORTED, (uint16_t)u, t->span.beg, p->last_end, m));
}

static void skip_line(parser_t* p) {
    while (!at_line_end(p)) advance(p);
}

// def name(params) [-> annotation]: suite
static void parse_def(parser_t* p) {
    const vial_token_t* t = cur(p);
    advance(p);
    const uint32_t m = mark(p);
    const vial_token_t* name = cur(p);
    if (name->kind == MD_VIAL_TOKEN_NAME) { advance(p); push(p, leaf(p, VIAL_N_NAME, 0, name)); }
    else { unexpected(p, "a function name"); push(p, error_node(p)); }

    const vial_token_t* open = cur(p);
    const uint32_t pm = mark(p);
    const uint32_t pbeg = open->span.beg;
    if (expect_op(p, VIAL_OP_LPAREN, "`(`")) {
        while (!is_op(p, VIAL_OP_RPAREN) && !at_line_end(p) && !p->panic) {
            const vial_token_t* pt = cur(p);
            if (tok_is_op(pt, VIAL_OP_STAR) || tok_is_op(pt, VIAL_OP_POW) || tok_is_op(pt, VIAL_OP_SLASH)) {
                refuse(p, "E0803", tok_span(pt), "%s", tok_is_op(pt, VIAL_OP_SLASH) ?
                    "positional-only parameters are not part of VIAL" : unsupported_message(VIAL_UNSUP_STARRED));
                advance(p);
                if (is_kind(p, MD_VIAL_TOKEN_NAME)) advance(p);
            } else if (pt->kind == MD_VIAL_TOKEN_NAME) {
                const uint32_t qm = mark(p);
                advance(p);
                push(p, leaf(p, VIAL_N_NAME, 0, pt));
                push(p, accept_op(p, VIAL_OP_COLON) ? parse_test(p) : 0);
                push(p, accept_op(p, VIAL_OP_ASSIGN) ? parse_test(p) : 0);
                push(p, make(p, VIAL_N_PARAM, 0, pt->span.beg, p->last_end, qm));
            } else {
                unexpected(p, "a parameter name");
                break;
            }
            if (!accept_op(p, VIAL_OP_COMMA)) break;
        }
        expect_close(p, open, VIAL_OP_RPAREN);
    }
    push(p, make(p, VIAL_N_PARAMS, 0, pbeg, p->last_end, pm));
    push(p, accept_op(p, VIAL_OP_ARROW) ? parse_test(p) : 0);
    push(p, parse_colon_suite(p, t));
    push(p, make(p, VIAL_N_DEF, 0, t->span.beg, p->last_end, m));
}

// if / elif chain, as nested IF nodes
static vial_node_id parse_if_chain(parser_t* p) {
    const vial_token_t* t = cur(p);
    advance(p);   // if or elif
    const uint32_t m = mark(p);
    push(p, parse_test(p));
    push(p, parse_colon_suite(p, t));
    if (is_kw(p, VIAL_KW_ELIF)) {
        p->panic = false;
        push(p, parse_if_chain(p));
    } else if (is_kw(p, VIAL_KW_ELSE)) {
        const vial_token_t* e = cur(p);
        advance(p);
        p->panic = false;
        push(p, parse_colon_suite(p, e));
    } else {
        push(p, 0);
    }
    return make(p, VIAL_N_IF, 0, t->span.beg, p->last_end, m);
}

static void parse_loop_else(parser_t* p) {
    if (!is_kw(p, VIAL_KW_ELSE)) return;
    const vial_token_t* e = cur(p);
    p->panic = false;
    refuse(p, "E0803", tok_span(e), "%s", unsupported_message(VIAL_UNSUP_LOOP_ELSE));
    advance(p);
    const uint32_t m = mark(p);
    push(p, parse_colon_suite(p, e));
    push(p, make(p, VIAL_N_UNSUPPORTED, VIAL_UNSUP_LOOP_ELSE, e->span.beg, p->last_end, m));
}

static void parse_for(parser_t* p) {
    const vial_token_t* t = cur(p);
    advance(p);
    const uint32_t m = mark(p);
    {
        const uint32_t tm = mark(p);
        const uint32_t tbeg = cur(p)->span.beg;
        push(p, parse_bitor(p));
        bool tuple = false;
        while (accept_op(p, VIAL_OP_COMMA)) {
            tuple = true;
            if (is_kw(p, VIAL_KW_IN)) break;
            push(p, parse_bitor(p));
        }
        if (tuple) push(p, make(p, VIAL_N_TUPLE, 0, tbeg, p->last_end, tm));
    }
    if (!accept_kw(p, VIAL_KW_IN)) { unexpected(p, "`in`"); push(p, error_node(p)); }
    else push(p, parse_testlist(p));
    push(p, parse_colon_suite(p, t));
    push(p, make(p, VIAL_N_FOR, 0, t->span.beg, p->last_end, m));
    parse_loop_else(p);
}

static void parse_while(parser_t* p) {
    const vial_token_t* t = cur(p);
    refuse(p, "E0801", tok_span(t), "`while` is not supported yet");
    advance(p);
    const uint32_t m = mark(p);
    push(p, parse_test(p));
    push(p, parse_colon_suite(p, t));
    push(p, make(p, VIAL_N_WHILE, 0, t->span.beg, p->last_end, m));
    parse_loop_else(p);
}

// Item assignment a[i] = x is not supported yet; reported on each target that is one
static void check_targets(parser_t* p, uint32_t first, uint32_t count) {
    for (uint32_t i = 0; i < count; ++i) {
        const vial_node_id id = p->stack[first + i];
        if (node_kind(p, id) == VIAL_N_SUBSCRIPT) {
            refuse(p, "E0801", p->tree.nodes[id].span, "%s", unsupported_message(VIAL_UNSUP_ITEM_ASSIGN));
        }
    }
}

static bool is_aug_assign(vial_op_t op) { return op >= VIAL_OP_PLUS_EQ && op <= VIAL_OP_RSHIFT_EQ; }

// One simple statement, without its line end
static void parse_small(parser_t* p) {
    const vial_token_t* t = cur(p);
    if (t->kind == MD_VIAL_TOKEN_KEYWORD) {
        switch (t->sub) {
        case VIAL_KW_PASS:     advance(p); push(p, leaf(p, VIAL_N_PASS, 0, t)); return;
        case VIAL_KW_BREAK:    advance(p); push(p, leaf(p, VIAL_N_BREAK, 0, t)); return;
        case VIAL_KW_CONTINUE: advance(p); push(p, leaf(p, VIAL_N_CONTINUE, 0, t)); return;
        case VIAL_KW_RETURN: {
            advance(p);
            const uint32_t m = mark(p);
            push(p, (at_line_end(p) || is_op(p, VIAL_OP_SEMICOLON)) ? 0 : parse_testlist(p));
            push(p, make(p, VIAL_N_RETURN, 0, t->span.beg, p->last_end, m));
            return;
        }
        case VIAL_KW_ASSERT: {
            advance(p);
            const uint32_t m = mark(p);
            push(p, parse_test(p));
            push(p, accept_op(p, VIAL_OP_COMMA) ? parse_test(p) : 0);
            push(p, make(p, VIAL_N_ASSERT, 0, t->span.beg, p->last_end, m));
            return;
        }
        case VIAL_KW_IMPORT: case VIAL_KW_FROM: case VIAL_KW_GLOBAL: case VIAL_KW_NONLOCAL: case VIAL_KW_DEL:
        case VIAL_KW_RAISE: {
            const vial_unsupported_t u = (t->sub == VIAL_KW_IMPORT || t->sub == VIAL_KW_FROM) ? VIAL_UNSUP_IMPORT :
                                         (t->sub == VIAL_KW_GLOBAL || t->sub == VIAL_KW_NONLOCAL) ? VIAL_UNSUP_GLOBAL :
                                         t->sub == VIAL_KW_DEL ? VIAL_UNSUP_DEL : VIAL_UNSUP_EXCEPTION;
            refuse(p, "E0803", tok_span(t), "%s", unsupported_message(u));
            while (!at_line_end(p) && !is_op(p, VIAL_OP_SEMICOLON)) advance(p);
            push(p, make(p, VIAL_N_UNSUPPORTED, (uint16_t)u, t->span.beg, p->last_end, mark(p)));
            return;
        }
        default: break;
        }
    }

    // An expression, an assignment or an augmented assignment
    const uint32_t m = mark(p);
    const vial_node_id first = parse_testlist(p);
    push(p, first);
    if (is_op(p, VIAL_OP_ASSIGN) && !p->panic) {
        while (accept_op(p, VIAL_OP_ASSIGN)) {
            push(p, parse_testlist(p));
            if (p->panic) break;
        }
        check_targets(p, m, mark(p) - m - 1);
        push(p, make(p, VIAL_N_ASSIGN, 0, t->span.beg, p->last_end, m));
        return;
    }
    const vial_token_t* op = cur(p);
    if (op->kind == MD_VIAL_TOKEN_OPERATOR && is_aug_assign((vial_op_t)op->sub) && !p->panic) {
        advance(p);
        push(p, parse_testlist(p));
        check_targets(p, m, 1);
        push(p, make(p, VIAL_N_AUG_ASSIGN, op->sub, t->span.beg, p->last_end, m));
        return;
    }
    if (tok_is_op(op, VIAL_OP_COLON) && !p->panic) {
        // x: int = 1
        refuse(p, "E0803", tok_span(op), "%s", unsupported_message(VIAL_UNSUP_ANNOTATION));
        advance(p);
        push(p, parse_test(p));
        if (accept_op(p, VIAL_OP_ASSIGN)) push(p, parse_testlist(p));
        push(p, make(p, VIAL_N_UNSUPPORTED, VIAL_UNSUP_ANNOTATION, t->span.beg, p->last_end, m));
        return;
    }
    push(p, make(p, VIAL_N_EXPR_STMT, 0, t->span.beg, p->last_end, m));
}

static void parse_stmt(parser_t* p) {
    p->panic = false;
    const vial_token_t* t = cur(p);

    if (t->kind == MD_VIAL_TOKEN_KEYWORD) {
        switch (t->sub) {
        case VIAL_KW_DEF:   parse_def(p); return;
        case VIAL_KW_IF:    push(p, parse_if_chain(p)); return;
        case VIAL_KW_FOR:   parse_for(p); return;
        case VIAL_KW_WHILE: parse_while(p); return;
        case VIAL_KW_CLASS: parse_refused_compound(p, VIAL_UNSUP_CLASS); return;
        case VIAL_KW_TRY:   parse_refused_compound(p, VIAL_UNSUP_EXCEPTION); return;
        case VIAL_KW_WITH:  parse_refused_compound(p, VIAL_UNSUP_WITH); return;
        case VIAL_KW_ASYNC:
            refuse(p, "E0803", tok_span(t), "%s", unsupported_message(VIAL_UNSUP_ASYNC));
            advance(p);
            parse_stmt(p);   // async def, async for: the rest as it is
            return;
        case VIAL_KW_ELIF: case VIAL_KW_ELSE: case VIAL_KW_EXCEPT: case VIAL_KW_FINALLY:
            unexpected(p, NULL);
            synchronize(p);
            return;
        default: break;
        }
    }
    if (tok_is_op(t, VIAL_OP_AT)) {
        // @decorator on its own line, then the definition as it is
        refuse(p, "E0803", tok_span(t), "%s", unsupported_message(VIAL_UNSUP_DECORATOR));
        skip_line(p);
        if (is_kind(p, MD_VIAL_TOKEN_NEWLINE)) advance(p);
        return;
    }

    // Simple statements, separated by ';'
    for (;;) {
        parse_small(p);
        if (p->panic) break;
        if (!accept_op(p, VIAL_OP_SEMICOLON)) break;
        if (at_line_end(p)) break;
    }
    if (p->panic) { synchronize(p); return; }
    if (is_kind(p, MD_VIAL_TOKEN_NEWLINE)) { advance(p); return; }
    if (is_kind(p, MD_VIAL_TOKEN_EOF) || is_kind(p, MD_VIAL_TOKEN_DEDENT)) return;
    unexpected(p, "end of line");
    synchronize(p);
}

vial_tree_t vial_parse(const vial_token_t* tokens, str_t src, vial_diags_t* diags, md_allocator_i* alloc) {
    parser_t p = {.src = src, .diags = diags, .alloc = alloc};
    for (size_t i = 0; i < md_array_size(tokens); ++i) {
        if (tokens[i].kind != MD_VIAL_TOKEN_COMMENT) md_array_push(p.toks, tokens[i], alloc);
    }
    if (md_array_size(p.toks) == 0 || md_array_last(p.toks)->kind != MD_VIAL_TOKEN_EOF) {
        vial_token_t eof = {.kind = MD_VIAL_TOKEN_EOF, .span = {(uint32_t)src.len, (uint32_t)src.len}};
        md_array_push(p.toks, eof, alloc);
    }
    p.num = md_array_size(p.toks);

    vial_node_t none = {0};
    md_array_push(p.tree.nodes, none, alloc);
    md_array_push(p.depth, 0, alloc);

    const uint32_t m = mark(&p);
    parse_stmts(&p, false);
    p.tree.root = make(&p, VIAL_N_MODULE, 0, 0, (uint32_t)src.len, m);
    vial_diags_sort(diags);

    md_array_free(p.toks, alloc);
    md_array_free(p.stack, alloc);
    md_array_free(p.depth, alloc);
    return p.tree;
}

// ### DUMP ###

const char* vial_node_kind_str(vial_node_kind_t kind) {
    static const char* names[VIAL_N_COUNT] = {
        [VIAL_N_NONE] = "-", [VIAL_N_MODULE] = "module", [VIAL_N_BLOCK] = "block", [VIAL_N_EXPR_STMT] = "expr",
        [VIAL_N_ASSIGN] = "assign", [VIAL_N_AUG_ASSIGN] = "augassign", [VIAL_N_DEF] = "def", [VIAL_N_PARAMS] = "params",
        [VIAL_N_PARAM] = "param", [VIAL_N_RETURN] = "return", [VIAL_N_IF] = "if", [VIAL_N_FOR] = "for",
        [VIAL_N_WHILE] = "while", [VIAL_N_BREAK] = "break", [VIAL_N_CONTINUE] = "continue", [VIAL_N_PASS] = "pass",
        [VIAL_N_ASSERT] = "assert", [VIAL_N_UNSUPPORTED] = "unsupported", [VIAL_N_ERROR] = "error",
        [VIAL_N_NAME] = "name", [VIAL_N_INT] = "int", [VIAL_N_FLOAT] = "float", [VIAL_N_STRING] = "string", [VIAL_N_QUANTITY] = "quantity", [VIAL_N_UNIT] = "unit",
        [VIAL_N_CONST] = "const", [VIAL_N_ELLIPSIS] = "ellipsis", [VIAL_N_UNARY] = "unary", [VIAL_N_NOT] = "not",
        [VIAL_N_BINARY] = "binary", [VIAL_N_BOOLOP] = "boolop", [VIAL_N_COMPARE] = "compare", [VIAL_N_CMPOP] = "cmpop",
        [VIAL_N_RANGE] = "range", [VIAL_N_TERNARY] = "ternary", [VIAL_N_CALL] = "call", [VIAL_N_KWARG] = "kwarg",
        [VIAL_N_ATTR] = "attr", [VIAL_N_SUBSCRIPT] = "subscript", [VIAL_N_SLICE] = "slice", [VIAL_N_TUPLE] = "tuple",
        [VIAL_N_LIST] = "list", [VIAL_N_LISTCOMP] = "listcomp", [VIAL_N_COMP_FOR] = "compfor",
    };
    return kind < VIAL_N_COUNT && names[kind] ? names[kind] : "?";
}

const char* vial_cmp_str(vial_cmp_t cmp) {
    switch (cmp) {
    case VIAL_CMP_LT: return "<";
    case VIAL_CMP_GT: return ">";
    case VIAL_CMP_LE: return "<=";
    case VIAL_CMP_GE: return ">=";
    case VIAL_CMP_EQ: return "==";
    case VIAL_CMP_NE: return "!=";
    case VIAL_CMP_IN: return "in";
    case VIAL_CMP_NOT_IN: return "not in";
    case VIAL_CMP_IS: return "is";
    case VIAL_CMP_IS_NOT: return "is not";
    default: return "?";
    }
}

static const char* unsupported_name(vial_unsupported_t u) {
    switch (u) {
    case VIAL_UNSUP_CLASS: return "class";        case VIAL_UNSUP_IMPORT: return "import";
    case VIAL_UNSUP_LAMBDA: return "lambda";      case VIAL_UNSUP_EXCEPTION: return "exception";
    case VIAL_UNSUP_GENERATOR: return "generator"; case VIAL_UNSUP_GLOBAL: return "global";
    case VIAL_UNSUP_DECORATOR: return "decorator"; case VIAL_UNSUP_WITH: return "with";
    case VIAL_UNSUP_ASYNC: return "async";        case VIAL_UNSUP_DICT: return "dict";
    case VIAL_UNSUP_SET: return "set";            case VIAL_UNSUP_STARRED: return "starred";
    case VIAL_UNSUP_FSTRING: return "fstring";    case VIAL_UNSUP_DEL: return "del";
    case VIAL_UNSUP_WALRUS: return "walrus";      case VIAL_UNSUP_LOOP_ELSE: return "loop-else";
    case VIAL_UNSUP_ANNOTATION: return "annotation"; case VIAL_UNSUP_BYTES: return "bytes";
    case VIAL_UNSUP_ITEM_ASSIGN: return "item-assign"; case VIAL_UNSUP_V1_RANGE: return "v1-range";
    default: return "?";
    }
}

static void dump_fmt(char** buf, md_allocator_i* alloc, const char* fmt, ...) {
    char tmp[512];
    va_list args;
    va_start(args, fmt);
    int len = vsnprintf(tmp, sizeof(tmp), fmt, args);
    va_end(args);
    if (len < 0) return;
    if ((size_t)len >= sizeof(tmp)) len = (int)sizeof(tmp) - 1;
    md_array_push_array(*buf, tmp, (size_t)len, alloc);
}

static void dump_node(char** buf, const vial_tree_t* t, vial_node_id id, str_t src, int depth, md_allocator_i* alloc) {
    for (int i = 0; i < depth; ++i) dump_fmt(buf, alloc, "  ");
    if (id == 0) { dump_fmt(buf, alloc, "-\n"); return; }
    const vial_node_t* n = &t->nodes[id];
    dump_fmt(buf, alloc, "%s", vial_node_kind_str((vial_node_kind_t)n->kind));
    switch (n->kind) {
    case VIAL_N_BINARY: case VIAL_N_UNARY: case VIAL_N_AUG_ASSIGN: dump_fmt(buf, alloc, " %s", vial_op_str((vial_op_t)n->sub)); break;
    case VIAL_N_BOOLOP: case VIAL_N_CONST: dump_fmt(buf, alloc, " %s", vial_kw_str((vial_kw_t)n->sub)); break;
    case VIAL_N_CMPOP: dump_fmt(buf, alloc, " %s", vial_cmp_str((vial_cmp_t)n->sub)); break;
    case VIAL_N_UNSUPPORTED: dump_fmt(buf, alloc, " %s", unsupported_name((vial_unsupported_t)n->sub)); break;
    case VIAL_N_NAME: case VIAL_N_INT: case VIAL_N_FLOAT: case VIAL_N_STRING: case VIAL_N_UNIT: {
        const uint32_t len = n->span.end - n->span.beg;
        dump_fmt(buf, alloc, " %.*s", (int)(len > 60 ? 60 : len), src.ptr + n->span.beg);
        break;
    }
    default: break;
    }
    dump_fmt(buf, alloc, " %u..%u\n", n->span.beg, n->span.end);
    for (uint32_t i = 0; i < n->child_count; ++i) {
        dump_node(buf, t, t->children[n->child_beg + i], src, depth + 1, alloc);
    }
}

str_t vial_tree_dump(const vial_tree_t* tree, str_t src, md_allocator_i* alloc) {
    char* buf = NULL;
    dump_node(&buf, tree, tree->root, src, 0, alloc);
    return (str_t){buf, md_array_size(buf)};
}

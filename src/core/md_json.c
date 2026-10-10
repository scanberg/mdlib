#include <core/md_json.h>

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_intrinsics.h>
#include <core/md_log.h>
#include <core/md_parse.h>

#include <string.h>

// The document is one flat array of tokens in document order, one per value and one per member key,
// with a member's key immediately in front of its value. Every token knows the index just past
// everything inside it (next), which is what makes stepping over a value O(1). Token 0 is a sentinel
// standing for NONE, so a lookup that finds nothing can still hand out an index.
//
// Parsing is one pass with no recursion and no stack of its own: while a container is open, its next
// field holds the index of the container around it, and is overwritten with the real next when the
// container closes. Nesting depth is therefore limited by nothing but the text.
//
// The tokens are allocated up front, sized by counting the characters every token but the root has
// in front of it - '[' or ',' for an element, '{' or ',' for a key, ':' for a member value - which
// is an upper bound (it counts those characters inside strings too) found by a loop the compiler
// vectorizes.

enum {
    TOK_LAST    = 1,    // the last element / member value of its container, or the root
    TOK_KEY     = 2,    // a member key
    TOK_ESCAPED = 4,    // a string with a backslash in it
    TOK_INT     = 8,    // a number written as an integer: no fraction, no exponent
};

typedef struct tok_t {
    uint32_t offset;    // where the value starts in the text; for a string, just after the quote
    uint32_t len;       // scalars: length of the text, strings without their quotes; containers: child count
    uint32_t next;      // index of the first token after this value and everything in it
    uint8_t  type;
    uint8_t  flags;
    uint16_t unused;
} tok_t;

struct md_json_t {
    const char* text;
    size_t      alloc_size;
    uint32_t    num_tok;
    uint32_t    unused;
    // tok_t tok[num_tok] follows
};

STATIC_ASSERT(sizeof(tok_t) == 16, "token is meant to be 16 bytes");
STATIC_ASSERT(sizeof(struct md_json_t) % 4 == 0, "tokens follow the header");

static inline const tok_t* doc_tok(const md_json_t* doc) {
    return (const tok_t*)(doc + 1);
}

static const tok_t none_tok = { 0, 0, 0, MD_JSON_TYPE_NONE, TOK_LAST, 0 };

static inline const tok_t* val_tok(md_json_val_t v) {
    if (!v.doc) return &none_tok;
    ASSERT(v.idx < v.doc->num_tok);
    return doc_tok(v.doc) + v.idx;
}

static inline md_json_val_t make_val(const md_json_t* doc, uint32_t idx) {
    md_json_val_t v = { doc, idx };
    return v;
}

static inline md_json_val_t none_val(void) {
    md_json_val_t v = { 0, 0 };
    return v;
}

// ---------------------------------------------------------------------------
// Tokenizer
// ---------------------------------------------------------------------------

// Character classes
enum {
    CC_WS    = 1,       // ' ' \t \n \r
    CC_DELIM = 2,       // what may follow a number or a literal: whitespace , ] }
    CC_STR   = 4,       // what ends a run of plain string bytes: " \ and control characters
    CC_DIGIT = 8,
    CC_HEX   = 16,
};

#define C_ CC_STR                       // control character
#define W_ (CC_WS | CC_DELIM | CC_STR)  // \t \n \r: whitespace, and control characters in a string
#define S_ (CC_WS | CC_DELIM)           // space
#define E_ CC_DELIM                     // , ] }
#define Q_ CC_STR                       // " and backslash
#define D_ (CC_DIGIT | CC_HEX)
#define H_ CC_HEX

static const uint8_t char_class[256] = {
    C_, C_, C_, C_, C_, C_, C_, C_, C_, W_, W_, C_, C_, W_, C_, C_,     // 0x00
    C_, C_, C_, C_, C_, C_, C_, C_, C_, C_, C_, C_, C_, C_, C_, C_,     // 0x10
    S_,  0, Q_,  0,  0,  0,  0,  0,  0,  0,  0,  0, E_,  0,  0,  0,     //  !"#$%&'()*+,-./
    D_, D_, D_, D_, D_, D_, D_, D_, D_, D_,  0,  0,  0,  0,  0,  0,     // 0123456789:;<=>?
     0, H_, H_, H_, H_, H_, H_,  0,  0,  0,  0,  0,  0,  0,  0,  0,     // @ABCDEFGHIJKLMNO
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, Q_, E_,  0,  0,     // PQRSTUVWXYZ[\]^_
     0, H_, H_, H_, H_, H_, H_,  0,  0,  0,  0,  0,  0,  0,  0,  0,     // `abcdefghijklmno
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, E_,  0,  0,     // pqrstuvwxyz{|}~
    // 0x80 - 0xFF: UTF-8 sequences, passed through
};

#undef C_
#undef W_
#undef S_
#undef E_
#undef Q_
#undef D_
#undef H_

#define CLASS(c) (char_class[(uint8_t)(c)])

// Pretty printed JSON is mostly indentation, so runs of spaces are taken 8 bytes at a time: the first
// byte that is not a space is the lowest nonzero byte of the word XOR 8 spaces (little endian).
static inline const char* skip_ws(const char* p, const char* end) {
    while (end - p >= 8) {
        uint64_t x;
        memcpy(&x, p, sizeof(x));
        x ^= 0x2020202020202020ull;
        if (x == 0) {
            p += 8;
            continue;
        }
        p += ctz64(x) >> 3;
        if (!(CLASS(*p) & CC_WS)) return p;
        ++p;
    }
    while (p < end && (CLASS(*p) & CC_WS)) ++p;
    return p;
}

static inline bool at_delim(const char* p, const char* end) {
    return p == end || (CLASS(*p) & CC_DELIM);
}

typedef struct scan_err_t {
    const char* msg;
    const char* at;
} scan_err_t;

// p is just after the opening quote. Gives the closing quote, or NULL with the error set.
static const char* scan_string(const char* p, const char* end, bool* escaped, scan_err_t* err) {
    for (;;) {
        while (p < end && !(CLASS(*p) & CC_STR)) ++p;
        if (p == end) {
            err->msg = "unterminated string";
            err->at  = p;
            return NULL;
        }
        const char c = *p;
        if (c == '"') {
            return p;
        }
        if (c != '\\') {
            err->msg = "control character in string, which must be escaped";
            err->at  = p;
            return NULL;
        }
        *escaped = true;
        const char* esc = p++;
        if (p == end) {
            err->msg = "unterminated string";
            err->at  = p;
            return NULL;
        }
        switch (*p) {
        case '"': case '\\': case '/': case 'b': case 'f': case 'n': case 'r': case 't':
            ++p;
            break;
        case 'u':
            if (end - p < 5 || !(CLASS(p[1]) & CLASS(p[2]) & CLASS(p[3]) & CLASS(p[4]) & CC_HEX)) {
                err->msg = "invalid \\u escape, which takes four hex digits";
                err->at  = esc;
                return NULL;
            }
            p += 5;
            break;
        default:
            err->msg = "invalid escape in string";
            err->at  = esc;
            return NULL;
        }
    }
}

static inline bool match_lit(const char* p, const char* end, const char* lit, size_t len) {
    return (size_t)(end - p) >= len && memcmp(p, lit, len) == 0;
}

// p is at '-' or a digit. Gives the end of the number, or NULL with the error set.
static const char* scan_number(const char* p, const char* end, bool* integer, scan_err_t* err) {
    const char* beg = p;
    *integer = true;
    if (*p == '-') {
        ++p;
        if (match_lit(p, end, "Infinity", 8)) {
            *integer = false;
            p += 8;
            goto delim;
        }
    }
    if (p == end) goto invalid;
    if (*p == '0') {
        ++p;
    } else if (CLASS(*p) & CC_DIGIT) {
        ++p;
        while (p < end && (CLASS(*p) & CC_DIGIT)) ++p;
    } else {
        goto invalid;
    }
    if (p < end && *p == '.') {
        *integer = false;
        ++p;
        if (p == end || !(CLASS(*p) & CC_DIGIT)) goto invalid;
        while (p < end && (CLASS(*p) & CC_DIGIT)) ++p;
    }
    if (p < end && (*p == 'e' || *p == 'E')) {
        *integer = false;
        ++p;
        if (p < end && (*p == '+' || *p == '-')) ++p;
        if (p == end || !(CLASS(*p) & CC_DIGIT)) goto invalid;
        while (p < end && (CLASS(*p) & CC_DIGIT)) ++p;
    }
delim:
    if (!at_delim(p, end)) goto invalid;
    return p;
invalid:
    err->msg = "invalid number";
    err->at  = beg;
    return NULL;
}

static inline unsigned is_separator(char c) {
    return (c == ',') | (c == ':') | (c == '[') | (c == '{');
}

static size_t token_bound(const char* p, const char* end) {
    // Counted into a byte per 64 byte block, which compilers vectorize at -O2 to 15-20 GB/s; a
    // size_t counter per byte is a quarter of that, as every compare has to be widened to it
    size_t count = 0;
    while (end - p >= 64) {
        uint8_t n = 0;
        for (int i = 0; i < 64; ++i) n += (uint8_t)is_separator(p[i]);
        count += n;
        p += 64;
    }
    for (; p < end; ++p) count += is_separator(*p);
    return count;
}

static void report_error(md_json_error_t* err, const char* msg, const char* beg, const char* at) {
    size_t line = 1;
    const char* line_beg = beg;
    for (const char* c = beg; c < at; ++c) {
        if (*c == '\n') {
            line += 1;
            line_beg = c + 1;
        }
    }
    const size_t offset = (size_t)(at - beg);
    const size_t column = (size_t)(at - line_beg) + 1;
    if (err) {
        err->message = msg;
        err->offset  = offset;
        err->line    = line;
        err->column  = column;
    } else {
        MD_LOG_ERROR("JSON: %s, at line %zu column %zu", msg, line, column);
    }
}

md_json_t* md_json_parse(str_t text, md_allocator_i* alloc, md_json_error_t* err) {
    ASSERT(alloc);

    const char* beg = text.ptr ? text.ptr : "";
    const size_t len = text.ptr ? text.len : 0;
    const char* end = beg + len;

    scan_err_t e = {0};
    md_json_t* doc = NULL;
    size_t alloc_size = 0;

    if (len >= UINT32_MAX) {
        e.msg = "text larger than 4 GiB";
        e.at  = beg;
        goto fail;
    }

    {
        const size_t max_tok = token_bound(beg, end) + 2;
        alloc_size = sizeof(md_json_t) + max_tok * sizeof(tok_t);
        doc = md_alloc(alloc, alloc_size);
        if (!doc) {
            e.msg = "out of memory";
            e.at  = beg;
            goto fail;
        }
        MEMSET(doc, 0, sizeof(md_json_t));
        tok_t* tok = (tok_t*)(doc + 1);
        tok[0] = none_tok;

        uint32_t n = 1;         // tokens so far
        uint32_t parent = 0;    // the innermost open container, 0 at the root
        uint32_t last = 0;      // the value completed most recently
        const char* p = beg;

        if (len >= 3 && (uint8_t)p[0] == 0xEF && (uint8_t)p[1] == 0xBB && (uint8_t)p[2] == 0xBF) {
            p += 3;
        }

value:
        p = skip_ws(p, end);
        if (p == end) {
            e.msg = "unexpected end of text, expected a value";
            e.at  = p;
            goto fail;
        }
        {
            ASSERT(n < max_tok);
            const uint32_t i = n++;
            tok_t* t = &tok[i];
            t->offset = (uint32_t)(p - beg);
            t->flags  = 0;
            t->unused = 0;
            const char c = *p;
            switch (c) {
            case '{':
            case '[':
                t->type = (c == '{') ? MD_JSON_TYPE_OBJECT : MD_JSON_TYPE_ARRAY;
                t->len  = 0;
                t->next = parent;   // the link to the enclosing container while this one is open
                parent  = i;
                p = skip_ws(p + 1, end);
                if (p < end && *p == (c == '{' ? '}' : ']')) {
                    ++p;
                    goto close;
                }
                if (c == '{') goto key;
                goto value;
            case '"': {
                bool escaped = false;
                const char* q = scan_string(p + 1, end, &escaped, &e);
                if (!q) goto fail;
                t->type   = MD_JSON_TYPE_STRING;
                t->offset = (uint32_t)(p + 1 - beg);
                t->len    = (uint32_t)(q - (p + 1));
                t->flags  = escaped ? TOK_ESCAPED : 0;
                p = q + 1;
                break;
            }
            case '-': case '0': case '1': case '2': case '3': case '4':
            case '5': case '6': case '7': case '8': case '9': {
                bool integer;
                const char* q = scan_number(p, end, &integer, &e);
                if (!q) goto fail;
                t->type  = MD_JSON_TYPE_NUMBER;
                t->len   = (uint32_t)(q - p);
                t->flags = integer ? TOK_INT : 0;
                p = q;
                break;
            }
            default: {
                static const struct {
                    const char* text;
                    uint8_t len;
                    uint8_t type;
                } lits[] = {
                    { "true",     4, MD_JSON_TYPE_BOOL   },
                    { "false",    5, MD_JSON_TYPE_BOOL   },
                    { "null",     4, MD_JSON_TYPE_NULL   },
                    { "NaN",      3, MD_JSON_TYPE_NUMBER },
                    { "Infinity", 8, MD_JSON_TYPE_NUMBER },
                };
                int k = 0;
                while (k < (int)ARRAY_SIZE(lits) && lits[k].text[0] != c) ++k;
                if (k == (int)ARRAY_SIZE(lits)) {
                    e.msg = "unexpected character, expected a value";
                    e.at  = p;
                    goto fail;
                }
                if (!match_lit(p, end, lits[k].text, lits[k].len) || !at_delim(p + lits[k].len, end)) {
                    e.msg = "invalid literal";
                    e.at  = p;
                    goto fail;
                }
                t->type = lits[k].type;
                t->len  = lits[k].len;
                p += lits[k].len;
                break;
            }
            }
            t->next = n;
            last = i;
        }

after:
        if (parent == 0) {
            p = skip_ws(p, end);
            if (p != end) {
                e.msg = "unexpected text after the root value";
                e.at  = p;
                goto fail;
            }
            goto done;
        }
        tok[parent].len += 1;
        p = skip_ws(p, end);
        {
            const bool in_object = tok[parent].type == MD_JSON_TYPE_OBJECT;
            if (p < end && *p == ',') {
                ++p;
                if (in_object) goto key;
                goto value;
            }
            if (p < end && *p == (in_object ? '}' : ']')) {
                ++p;
                goto close;
            }
            e.msg = (p == end) ? (in_object ? "unexpected end of text, expected ',' or '}'" : "unexpected end of text, expected ',' or ']'")
                               : (in_object ? "expected ',' or '}'" : "expected ',' or ']'");
            e.at  = p;
            goto fail;
        }

close:
        {
            const uint32_t c = parent;
            parent = tok[c].next;
            tok[c].next = n;
            if (tok[c].len) tok[last].flags |= TOK_LAST;
            last = c;
            goto after;
        }

key:
        p = skip_ws(p, end);
        if (p == end || *p != '"') {
            e.msg = (p == end) ? "unexpected end of text, expected a key" : "expected a string key";
            e.at  = p;
            goto fail;
        }
        {
            ASSERT(n < max_tok);
            const uint32_t i = n++;
            tok_t* t = &tok[i];
            bool escaped = false;
            const char* q = scan_string(p + 1, end, &escaped, &e);
            if (!q) goto fail;
            t->type   = MD_JSON_TYPE_STRING;
            t->flags  = TOK_KEY | (escaped ? TOK_ESCAPED : 0);
            t->unused = 0;
            t->offset = (uint32_t)(p + 1 - beg);
            t->len    = (uint32_t)(q - (p + 1));
            t->next   = n;
            p = skip_ws(q + 1, end);
        }
        if (p == end || *p != ':') {
            e.msg = (p == end) ? "unexpected end of text, expected ':'" : "expected ':' after a key";
            e.at  = p;
            goto fail;
        }
        ++p;
        goto value;

done:
        tok[1].flags |= TOK_LAST;
        doc->text       = beg;
        doc->alloc_size = alloc_size;
        doc->num_tok    = n;
        return doc;
    }

fail:
    if (doc) md_free(alloc, doc, alloc_size);
    report_error(err, e.msg, beg, e.at);
    return NULL;
}

void md_json_free(md_json_t* doc, md_allocator_i* alloc) {
    ASSERT(alloc);
    if (doc) md_free(alloc, doc, doc->alloc_size);
}

// ---------------------------------------------------------------------------
// Navigation
// ---------------------------------------------------------------------------

md_json_val_t md_json_root(const md_json_t* doc) {
    return doc ? make_val(doc, 1) : none_val();
}

md_json_type_t md_json_type(md_json_val_t v) {
    return (md_json_type_t)val_tok(v)->type;
}

size_t md_json_count(md_json_val_t v) {
    const tok_t* t = val_tok(v);
    return (t->type == MD_JSON_TYPE_ARRAY || t->type == MD_JSON_TYPE_OBJECT) ? t->len : 0;
}

md_json_val_t md_json_first(md_json_val_t v) {
    const tok_t* t = val_tok(v);
    if (t->type == MD_JSON_TYPE_ARRAY  && t->len) return make_val(v.doc, v.idx + 1);
    if (t->type == MD_JSON_TYPE_OBJECT && t->len) return make_val(v.doc, v.idx + 2);
    return none_val();
}

md_json_val_t md_json_next(md_json_val_t v) {
    const tok_t* t = val_tok(v);
    if (t->flags & (TOK_LAST | TOK_KEY)) return none_val();
    uint32_t i = t->next;
    if (doc_tok(v.doc)[i].flags & TOK_KEY) i += 1;
    return make_val(v.doc, i);
}

md_json_val_t md_json_at(md_json_val_t v, size_t i) {
    const tok_t* t = val_tok(v);
    if (t->type == MD_JSON_TYPE_ARRAY && i < t->len) {
        const tok_t* tok = doc_tok(v.doc);
        uint32_t k = v.idx + 1;
        for (size_t n = 0; n < i; ++n) k = tok[k].next;
        return make_val(v.doc, k);
    }
    if (t->type == MD_JSON_TYPE_OBJECT && i < t->len) {
        const tok_t* tok = doc_tok(v.doc);
        uint32_t k = v.idx + 1;     // the key
        for (size_t n = 0; n < i; ++n) k = tok[k + 1].next;
        return make_val(v.doc, k + 1);
    }
    return none_val();
}

md_json_val_t md_json_key(md_json_val_t v) {
    const tok_t* t = val_tok(v);
    // The token in front of a member value is its key, and nothing else has a key in front of it
    if (t->type == MD_JSON_TYPE_NONE || (t->flags & TOK_KEY) || v.idx < 2) return none_val();
    if (doc_tok(v.doc)[v.idx - 1].flags & TOK_KEY) return make_val(v.doc, v.idx - 1);
    return none_val();
}

// ---------------------------------------------------------------------------
// Strings
// ---------------------------------------------------------------------------

static inline str_t tok_text(const md_json_t* doc, const tok_t* t) {
    str_t s = { doc->text + t->offset, t->len };
    return s;
}

static inline uint32_t hex4(const char* p) {
    uint32_t v = 0;
    for (int i = 0; i < 4; ++i) {
        const char c = p[i];
        const uint32_t d = (c <= '9') ? (uint32_t)(c - '0') : (uint32_t)((c | 0x20) - 'a' + 10);
        v = (v << 4) | d;
    }
    return v;
}

static inline size_t utf8_encode(char out[4], uint32_t cp) {
    if (cp < 0x80) {
        out[0] = (char)cp;
        return 1;
    }
    if (cp < 0x800) {
        out[0] = (char)(0xC0 | (cp >> 6));
        out[1] = (char)(0x80 | (cp & 0x3F));
        return 2;
    }
    if (cp < 0x10000) {
        out[0] = (char)(0xE0 | (cp >> 12));
        out[1] = (char)(0x80 | ((cp >> 6) & 0x3F));
        out[2] = (char)(0x80 | (cp & 0x3F));
        return 3;
    }
    out[0] = (char)(0xF0 | (cp >> 18));
    out[1] = (char)(0x80 | ((cp >> 12) & 0x3F));
    out[2] = (char)(0x80 | ((cp >> 6) & 0x3F));
    out[3] = (char)(0x80 | (cp & 0x3F));
    return 4;
}

// Decodes one character of a validated string body at *p into out, advancing *p. A character is an
// escape, or a UTF-8 lead byte with the continuation bytes after it, so a copy cut between two of
// them never splits a sequence.
static inline size_t decode_char(char out[4], const char** p, const char* end) {
    const char* c = *p;
    if (*c != '\\') {
        const uint8_t lead = (uint8_t)*c;
        size_t n = 1;
        if (lead >= 0xC0) {
            const size_t want = lead >= 0xF0 ? 4 : lead >= 0xE0 ? 3 : 2;
            while (n < want && c + n < end && ((uint8_t)c[n] & 0xC0) == 0x80) ++n;
        }
        memcpy(out, c, n);
        *p = c + n;
        return n;
    }
    const char e = c[1];
    *p = c + 2;
    switch (e) {
    case 'b': out[0] = '\b'; return 1;
    case 'f': out[0] = '\f'; return 1;
    case 'n': out[0] = '\n'; return 1;
    case 'r': out[0] = '\r'; return 1;
    case 't': out[0] = '\t'; return 1;
    case 'u': break;
    default:  out[0] = e;    return 1;     // " \ /
    }
    uint32_t cp = hex4(c + 2);
    *p = c + 6;
    if (0xD800 <= cp && cp <= 0xDBFF) {
        // A high surrogate takes the low one after it, when there is one
        if (end - *p >= 6 && (*p)[0] == '\\' && (*p)[1] == 'u') {
            const uint32_t lo = hex4(*p + 2);
            if (0xDC00 <= lo && lo <= 0xDFFF) {
                cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
                *p += 6;
                return utf8_encode(out, cp);
            }
        }
        cp = 0xFFFD;
    } else if (0xDC00 <= cp && cp <= 0xDFFF) {
        cp = 0xFFFD;
    }
    return utf8_encode(out, cp);
}

static bool string_eq(const md_json_t* doc, const tok_t* t, str_t str) {
    const str_t raw = tok_text(doc, t);
    if (!(t->flags & TOK_ESCAPED)) {
        return raw.len == str.len && (raw.len == 0 || memcmp(raw.ptr, str.ptr, raw.len) == 0);
    }
    const char* p   = raw.ptr;
    const char* end = raw.ptr + raw.len;
    size_t pos = 0;
    while (p < end) {
        char buf[4];
        const size_t n = decode_char(buf, &p, end);
        if (pos + n > str.len || memcmp(buf, str.ptr + pos, n) != 0) return false;
        pos += n;
    }
    return pos == str.len;
}

md_json_val_t md_json_get(md_json_val_t v, str_t key) {
    const tok_t* t = val_tok(v);
    if (t->type != MD_JSON_TYPE_OBJECT) return none_val();
    const tok_t* tok = doc_tok(v.doc);
    uint32_t k = v.idx + 1;
    uint32_t found = 0;
    for (uint32_t m = 0; m < t->len; ++m) {
        if (string_eq(v.doc, &tok[k], key)) found = k + 1;
        k = tok[k + 1].next;
    }
    return found ? make_val(v.doc, found) : none_val();
}

str_t md_json_string_raw(md_json_val_t v) {
    const tok_t* t = val_tok(v);
    if (t->type != MD_JSON_TYPE_STRING) {
        str_t empty = { 0, 0 };
        return empty;
    }
    return tok_text(v.doc, t);
}

bool md_json_string_eq(md_json_val_t v, str_t str) {
    const tok_t* t = val_tok(v);
    return t->type == MD_JSON_TYPE_STRING && string_eq(v.doc, t, str);
}

size_t md_json_string_copy(char* buf, size_t cap, md_json_val_t v) {
    ASSERT(buf || cap == 0);
    if (cap == 0) return 0;
    const tok_t* t = val_tok(v);
    size_t len = 0;
    if (t->type == MD_JSON_TYPE_STRING) {
        const str_t raw = tok_text(v.doc, t);
        const char* p   = raw.ptr;
        const char* end = raw.ptr + raw.len;
        while (p < end) {
            char ch[4];
            const size_t n = decode_char(ch, &p, end);
            if (len + n > cap - 1) break;
            memcpy(buf + len, ch, n);
            len += n;
        }
    }
    buf[len] = '\0';
    return len;
}

str_t md_json_string(md_json_val_t v, md_allocator_i* alloc) {
    ASSERT(alloc);
    str_t result = { 0, 0 };
    const tok_t* t = val_tok(v);
    if (t->type != MD_JSON_TYPE_STRING) return result;
    const str_t raw = tok_text(v.doc, t);
    if (!(t->flags & TOK_ESCAPED)) return str_copy(raw, alloc);

    // Escapes only ever shrink, so the raw length is enough
    char* buf = md_alloc(alloc, raw.len + 1);
    if (!buf) return result;
    const char* p   = raw.ptr;
    const char* end = raw.ptr + raw.len;
    size_t len = 0;
    while (p < end) {
        len += decode_char(buf + len, &p, end);
    }
    buf[len] = '\0';
    result.ptr = buf;
    result.len = len;
    return result;
}

// ---------------------------------------------------------------------------
// Numbers
// ---------------------------------------------------------------------------

bool md_json_f64(double* out, md_json_val_t v) {
    ASSERT(out);
    const tok_t* t = val_tok(v);
    if (t->type != MD_JSON_TYPE_NUMBER) return false;
    double value = 0.0;
    const size_t n = md_parse_f64(&value, tok_text(v.doc, t));
    ASSERT(n == t->len);
    (void)n;
    *out = value;
    return true;
}

bool md_json_i64(int64_t* out, md_json_val_t v) {
    ASSERT(out);
    const tok_t* t = val_tok(v);
    if (t->type != MD_JSON_TYPE_NUMBER) return false;
    const str_t text = tok_text(v.doc, t);

    if (t->flags & TOK_INT) {
        int64_t value = 0;
        md_parse_i64(&value, text);
        // md_parse_i64 saturates; the two limits themselves are the only literals that may give them
        if (value == INT64_MAX && !str_eq(text, STR_LIT("9223372036854775807")))  return false;
        if (value == INT64_MIN && !str_eq(text, STR_LIT("-9223372036854775808"))) return false;
        *out = value;
        return true;
    }

    double d = 0.0;
    md_parse_f64(&d, text);
    // On the bits, as this file is built with fast math: finite and |d| < 2^53, where every double
    // is a distinct integer or has a fraction
    uint64_t bits;
    memcpy(&bits, &d, sizeof(bits));
    const uint64_t mag = bits & 0x7FFFFFFFFFFFFFFFull;
    if (mag >= 0x4340000000000000ull) return false;     // 2^53
    const int64_t i = (int64_t)d;
    if ((double)i != d) return false;
    *out = i;
    return true;
}

bool md_json_bool(bool* out, md_json_val_t v) {
    ASSERT(out);
    const tok_t* t = val_tok(v);
    if (t->type != MD_JSON_TYPE_BOOL) return false;
    *out = v.doc->text[t->offset] == 't';
    return true;
}

size_t md_json_extract_f64(double* out, size_t cap, md_json_val_t v) {
    ASSERT(out || cap == 0);
    const tok_t* t = val_tok(v);
    if (t->type != MD_JSON_TYPE_ARRAY) return 0;
    const size_t count = MIN((size_t)t->len, cap);
    const tok_t* tok = doc_tok(v.doc);
    uint32_t k = v.idx + 1;
    size_t i = 0;
    for (; i < count; ++i) {
        const tok_t* e = &tok[k];
        if (e->type != MD_JSON_TYPE_NUMBER) break;
        md_parse_f64(&out[i], tok_text(v.doc, e));
        k = e->next;
    }
    return i;
}

size_t md_json_extract_f32(float* out, size_t cap, md_json_val_t v) {
    ASSERT(out || cap == 0);
    const tok_t* t = val_tok(v);
    if (t->type != MD_JSON_TYPE_ARRAY) return 0;
    const size_t count = MIN((size_t)t->len, cap);
    const tok_t* tok = doc_tok(v.doc);
    uint32_t k = v.idx + 1;
    size_t i = 0;
    for (; i < count; ++i) {
        const tok_t* e = &tok[k];
        if (e->type != MD_JSON_TYPE_NUMBER) break;
        double d = 0.0;
        md_parse_f64(&d, tok_text(v.doc, e));
        out[i] = (float)d;
        k = e->next;
    }
    return i;
}
